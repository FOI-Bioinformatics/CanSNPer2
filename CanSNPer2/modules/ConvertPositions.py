'''
Conversion of SNP positions from the format before 2.3.0 to the current format.

Before 2.3.0 a SNP was stored with one position along the sequences of its reference fasta
file joined in file order (the order NCBI uses). The current format stores the sequence (the
first word of its fasta header) and the 1-based position on that sequence. A joined position
p falls on the first sequence whose end is at or beyond p, and its position there is p minus
the lengths of the sequences before it.

The conversion depends on the order of the sequences in the reference file, so every converted
SNP is checked: the base at its new position must be its derived or its ancestral base. A row
that fails this check, or whose position is beyond the end of the joined sequences, is reported
with its location, and nothing is written when any row is refused.

Two kinds of input are converted:
- a tab separated source file with genome, position, derived_base and ancestral_base columns
  (the _snp.txt and _canSNPs.txt files); the position column is replaced by sequence and
  sequence_position, comment lines, blank lines and the other columns are copied unchanged;
- a database, whose snp_annotation table is rebuilt with the sequence and sequence_position
  columns in a copy of the database.
'''

import bisect
import io
import os
import shutil
import sqlite3
import tempfile
import logging

from CanSNPer2.modules.ReferenceCheck import read_sequences, listing, FastaFileError
from CanSNPer2.modules.DatabaseConnection import check_database, ConnectionError, CANSNP_SOURCE_COLUMNS
from CanSNPer2.modules.InitiateDatabase import snp_annotation_sql
from CanSNPer2.modules.FileMode import set_default_mode
from CanSNPer2.modules.SchemeInfo import write_scheme_info

logger = logging.getLogger(__name__)

'''refused rows listed in the message of a ConversionError before it is cut short'''
REPORTED = 50

'''columns of a source file in the format before 2.3.0 that the conversion reads'''
REQUIRED_COLUMNS = ("genome", "position", "derived_base", "ancestral_base")


class ConversionError(Exception):
	'''an input that cannot be converted; nothing has been written'''
	pass


def read_text(path, newline=None):
	'''the text of a SNP or references file as a file object. A UTF-8 byte order mark is
		removed; a gzip-compressed file or one that is not UTF-8 is refused with the name of
		the file, where the decoder named only the byte it stopped at'''
	with open(path, "rb") as fin:
		data = fin.read()
	if data[:2] == b"\x1f\x8b":
		raise ConversionError("{path} is gzip-compressed, decompress it first".format(path=path))
	try:
		return io.StringIO(data.decode("utf-8-sig"), newline=newline)
	except UnicodeDecodeError as e:
		raise ConversionError("{path} is not UTF-8 text ({error}), save it as UTF-8 and run again".format(path=path, error=e))


class JoinedReference(object):
	'''the sequences of a reference fasta file in file order, addressed by joined position'''
	def __init__(self, path):
		self.path = path
		try:
			self.records = read_sequences(path)
		except FastaFileError as e:
			raise ConversionError(str(e))
		if not self.records:
			raise ConversionError("The reference file {path} contains no sequences".format(path=path))
		self.starts = []  ## 0-based joined start of each sequence
		total = 0
		for _name, sequence in self.records:
			self.starts.append(total)
			total += len(sequence)
		self.length = total

	def locate(self, position):
		'''(sequence name, position on that sequence, base) of a 1-based joined position, or
			None when the position is beyond the end of the joined sequences'''
		if position < 1 or position > self.length:
			return None
		i = bisect.bisect_right(self.starts, position - 1) - 1
		name, sequence = self.records[i]
		local = position - self.starts[i]
		return name, local, sequence[local - 1]


def check_target(source, target):
	'''a conversion never writes over its input or over an existing file'''
	if os.path.exists(target):
		if os.path.realpath(source) == os.path.realpath(target) or (os.path.exists(source) and os.path.samefile(source, target)):
			raise ConversionError("The output {target} is the input, give another output folder with -o".format(target=target))
		raise ConversionError("The output {target} already exists, remove it or give another output folder with -o".format(target=target))


def _temporary(target):
	'''a new, empty file next to target that the output is written to before it is renamed to
		target; mkstemp never opens an existing file, so no other file is overwritten'''
	directory = os.path.dirname(target)
	if directory:
		os.makedirs(directory, exist_ok=True)
	handle, path = tempfile.mkstemp(prefix=".{name}.".format(name=os.path.basename(target)),
		suffix=".part", dir=directory or ".")
	os.close(handle)
	set_default_mode(path) ## the output gets the mode of open(), not the 0600 of mkstemp
	return path


def _remove(path):
	if os.path.exists(path):
		os.remove(path)


class _References(object):
	'''the JoinedReference of each genome, read once'''
	def __init__(self, genome_files):
		self.genome_files = genome_files
		self.loaded = {}

	def get(self, genome):
		if genome not in self.loaded:
			path = self.genome_files.get(genome)
			if path is None:
				raise KeyError(genome)
			self.loaded[genome] = JoinedReference(path)
		return self.loaded[genome]


def _convert_position(references, genome, position, derived, ancestral):
	'''((sequence, sequence_position), None) for a valid SNP, otherwise (None, reason)'''
	try:
		reference = references.get(genome)
	except KeyError:
		return None, "no reference file is given for genome {genome}".format(genome=genome)
	try:
		joined = int(position)
	except (TypeError, ValueError):
		return None, "the position {position} is not a whole number".format(position=position)
	if joined < 1:
		return None, "the joined position {pos} is not a positive position, positions count from 1".format(pos=joined)
	located = reference.locate(joined)
	if located is None:
		return None, ("joined position {pos} is beyond the end of the joined sequences of {path} "
			"({length} bases)").format(pos=joined, path=reference.path, length=reference.length)
	sequence, local, found = located
	if found not in (str(derived).upper(), str(ancestral).upper()):
		return None, ("joined position {pos} is {sequence}:{local}, where {path} has {found}, neither the derived "
			"base {derived} nor the ancestral base {ancestral}").format(pos=joined, sequence=sequence, local=local,
			path=reference.path, found=found, derived=derived, ancestral=ancestral)
	return (sequence, local), None


def _refusal(what, refused):
	shown = refused[:REPORTED]
	more = "" if len(refused) <= REPORTED else "\n  ... and {n} more".format(n=len(refused) - REPORTED)
	return ("{what} was not converted, {n} SNP(s) do not match the reference. The reference sequences must be in "
		"the order the joined positions were counted in (the order of the NCBI assembly file). Nothing was written.\n  "
		"{rows}{more}").format(what=what, n=len(refused), rows="\n  ".join(shown), more=more)


def _split_ending(line):
	if line.endswith("\r\n"):
		return line[:-2], "\r\n"
	if line.endswith("\n"):
		return line[:-1], "\n"
	return line, ""


def convert_table(source, target, genome_files, sequences=None):
	'''Convert a tab separated SNP file from joined positions to sequence and sequence_position
		and write it to target. genome_files: {genome: reference fasta path}. Returns the number
		of converted SNPs; the names of the sequences the SNPs are placed on are added to the
		set sequences when it is given. Raises ConversionError, and writes nothing, when any row
		is refused.'''
	check_target(source, target)
	references = _References(genome_files)
	output = []
	refused = []
	placed = set()
	header = None
	count = 0
	with read_text(source, newline="") as fin:
		for number, raw in enumerate(fin, start=1):
			line, ending = _split_ending(raw)
			if line.strip() == "" or line.startswith("#"):
				output.append(raw)
				continue
			fields = line.split("\t")
			if header is None:
				header = fields
				missing = [column for column in REQUIRED_COLUMNS if column not in header]
				if "sequence_position" in header:
					raise ConversionError("{source} already has a sequence_position column, it is in the current format".format(source=source))
				if missing:
					raise ConversionError("{source} has no column(s) {missing}, a header row with {required} is required".format(
						source=source, missing=", ".join(missing), required=", ".join(REQUIRED_COLUMNS)))
				column = {name: i for i, name in enumerate(header)}
				pos_i = column["position"]
				name_i = column.get("snp_id", column.get("node"))
				output.append("\t".join(header[:pos_i] + ["sequence", "sequence_position"] + header[pos_i + 1:]) + ending)
				continue
			if len(fields) < len(header):
				refused.append("{source} line {n}: the row has {f} columns, the header {h}".format(
					source=source, n=number, f=len(fields), h=len(header)))
				continue
			name = fields[name_i] if name_i is not None else "row"
			converted, reason = _convert_position(references, fields[column["genome"]], fields[pos_i],
				fields[column["derived_base"]], fields[column["ancestral_base"]])
			if converted is None:
				refused.append("{source} line {n}: {name}: {reason}".format(source=source, n=number, name=name, reason=reason))
				continue
			placed.add(converted[0])
			count += 1
			for key in ("derived_base", "ancestral_base"):
				fields[column[key]] = fields[column[key]].upper() ## --create and --import refuse lower case bases
			output.append("\t".join(fields[:pos_i] + [converted[0], str(converted[1])] + fields[pos_i + 1:]) + ending)
	if header is None:
		raise ConversionError("{source} has no header row".format(source=source))
	if refused:
		raise ConversionError(_refusal(source, refused))
	partial = _temporary(target)
	try:
		with open(partial, "w", newline="") as out:
			out.write("".join(output))
		os.replace(partial, target)
	finally:
		_remove(partial)
	if sequences is not None:
		sequences.update(placed)
	return count


def canSNPs_target(outdir, source):
	'''the path a converted SNP source file is written to: the name of the source with a
		trailing _snp or _canSNPs replaced by _canSNPs'''
	name = os.path.basename(source)
	stem = name[:-len(".txt")] if name.endswith(".txt") else name
	for ending in ("_snp", "_canSNPs"):
		if stem.endswith(ending):
			stem = stem[:-len(ending)]
			break
	return os.path.join(outdir, stem + "_canSNPs.txt")


def upgrade_table(source, target):
	'''Rewrite a SNP file of 2.3.0, or a canSNP file without the rank, score and typed columns,
		as a canSNP file. snp_id becomes node, the strain column is dropped, and the columns the
		file lacks are added: rank 1 for a SNP file, an empty score, and typed yes for the canSNP
		of rank 1 of each node. A canSNP file without a reference column is given NA, which is
		what the published schemes record for a canSNP whose source is not known.
		Comment and blank lines are copied, and any trailing columns of the source are kept after
		the typed column. Returns the number of canSNPs; raises ConversionError, and writes
		nothing, when a row cannot be read.'''
	check_target(source, target)
	header = None
	output = []
	refused = []
	count = 0
	with read_text(source, newline="") as fin:
		for number, raw in enumerate(fin, start=1):
			line, ending = _split_ending(raw)
			if line.strip() == "" or line.startswith("#"):
				output.append(raw)
				continue
			fields = line.split("\t")
			if header is None:
				header = fields
				if "sequence_position" not in header:
					raise ConversionError("{source} has no sequence_position column, convert its positions "
						"first".format(source=source))
				if "typed" in header:
					raise ConversionError("{source} already has a typed column, it is in the current format, "
						"nothing to upgrade".format(source=source))
				if "node" not in header and "snp_id" not in header:
					raise ConversionError("{source} has no node column, a header row with {required} is "
						"required".format(source=source, required=", ".join(CANSNP_SOURCE_COLUMNS)))
				column = {name: i for i, name in enumerate(header)}
				missing = [name for name in ("genome", "sequence", "sequence_position", "derived_base",
					"ancestral_base") if name not in column]
				if missing:
					raise ConversionError("{source} has no column(s) {missing}".format(
						source=source, missing=", ".join(missing)))
				kept = [name for name in header if name not in CANSNP_SOURCE_COLUMNS
					and name not in ("snp_id", "strain")]
				output.append("\t".join(CANSNP_SOURCE_COLUMNS + tuple(kept)) + ending)
				continue
			if len(fields) < len(header):
				refused.append("{source} line {n}: the row has {f} columns, the header {h}".format(
					source=source, n=number, f=len(fields), h=len(header)))
				continue
			value = lambda name, default="": fields[column[name]] if name in column else default
			node = value("node") or value("snp_id")
			rank = value("rank", "1") or "1"
			row = [node, rank, value("reference", "NA") or "NA", value("genome"), value("sequence"),
				value("sequence_position"), value("derived_base").upper(), value("ancestral_base").upper(),
				value("score"), "yes" if rank == "1" else "no"]
			output.append("\t".join(row + [value(name) for name in kept]) + ending)
			count += 1
	if header is None:
		raise ConversionError("{source} has no header row".format(source=source))
	if refused:
		raise ConversionError(_refusal(source, refused))
	partial = _temporary(target)
	try:
		with open(partial, "w", newline="") as out:
			out.write("".join(output))
		os.replace(partial, target)
	finally:
		_remove(partial)
	return count


def upgrade_database(source_db, target_db, scheme_info=None):
	'''Copy a database of 2.3.0, which holds one canSNP per node, to target_db and rebuild its
		snp_annotation table in the current format: every canSNP becomes rank 1 of its node,
		typed, without a score. No reference file is needed, no position moves.
		Returns the number of canSNPs; raises ConversionError and leaves no target_db.'''
	try:
		check_database(source_db, current_format=False)
	except ConnectionError as e:
		raise ConversionError(str(e))
	check_target(source_db, target_db)
	if database_stage(source_db) == "joined":
		raise ConversionError("The database {db} stores positions along the joined reference sequences; convert "
			"those first, which needs the reference files".format(db=source_db))
	try:
		rows = sqlite3.connect(source_db).execute('''SELECT node_id, snp_id, sequence, sequence_position,
			ancestral_base, derived_base, reference, date, genome_i FROM snp_annotation ORDER BY node_id''').fetchall()
	except sqlite3.Error as e:
		raise ConversionError("The canSNPs of {db} could not be read: {error}".format(db=source_db, error=e))
	partial = _temporary(target_db)
	try:
		shutil.copyfile(source_db, partial)
		conn = sqlite3.connect(partial)
		try:
			conn.execute(snp_annotation_sql("snp_annotation_upgraded"))
			conn.executemany('''INSERT INTO snp_annotation_upgraded (node_id, snp_id, sequence, sequence_position,
				ancestral_base, derived_base, reference, date, genome_i, rank, score, typed)
				VALUES (?,?,?,?,?,?,?,?,?,1,NULL,1)''', rows)
			conn.execute("DROP TABLE snp_annotation")
			conn.execute("ALTER TABLE snp_annotation_upgraded RENAME TO snp_annotation")
			conn.commit()
			if scheme_info:
				write_scheme_info(conn, scheme_info)
			conn.execute("VACUUM")
		finally:
			conn.close()
		try:
			check_database(partial)
		except ConnectionError as e:
			raise ConversionError("The upgraded copy of {db} failed the database check, nothing was written: "
				"{error}".format(db=source_db, error=e))
		os.replace(partial, target_db)
	except sqlite3.Error as e:
		raise ConversionError("The database {db} could not be upgraded, nothing was written: {error}".format(
			db=source_db, error=e))
	finally:
		_remove(partial)
	return len(rows)


def convert_database(source_db, target_db, genome_files, sequences=None, scheme_info=None):
	'''Copy a database from the format before 2.3.0 to target_db and rebuild its snp_annotation
		table with sequence and sequence_position columns (create, copy, drop and rename; ALTER
		TABLE DROP COLUMN needs sqlite 3.35). genome_files: {genome: reference fasta path}.
		scheme_info: keys written to the scheme_info table of the converted copy.
		Returns the number of converted SNPs. Raises ConversionError, and leaves no target_db,
		when any SNP is refused.'''
	try:
		check_database(source_db, current_format=False)
	except ConnectionError as e:
		raise ConversionError(str(e))
	check_target(source_db, target_db)
	try:
		rows = _joined_rows(source_db)
	except sqlite3.Error as e:
		raise ConversionError("The SNPs of {db} could not be read: {error}".format(db=source_db, error=e))
	references = _References(genome_files)
	refused = []
	converted_rows = []
	placed = set()
	for node_id, snp_id, position, ancestral, derived, reference, date, genome_i, genome in rows:
		if genome is None:
			refused.append("{db}: {snp} (node {node}): its genome_i {gi} is not in snp_references".format(
				db=source_db, snp=snp_id, node=node_id, gi=genome_i))
			continue
		converted, reason = _convert_position(references, genome, position, derived, ancestral)
		if converted is None:
			refused.append("{db}: {snp} (node {node}): {reason}".format(db=source_db, snp=snp_id, node=node_id, reason=reason))
			continue
		placed.add(converted[0])
		converted_rows.append((node_id, snp_id, converted[0], converted[1], str(ancestral).upper(), str(derived).upper(), reference, date, genome_i))
	if refused:
		raise ConversionError(_refusal(source_db, refused))
	partial = _temporary(target_db)
	try:
		shutil.copyfile(source_db, partial)
		conn = sqlite3.connect(partial)
		try:
			conn.execute(snp_annotation_sql("snp_annotation_converted"))
			conn.executemany('''INSERT INTO snp_annotation_converted (node_id, snp_id, sequence, sequence_position,
				ancestral_base, derived_base, reference, date, genome_i) VALUES (?,?,?,?,?,?,?,?,?)''', converted_rows)
			conn.execute("DROP TABLE snp_annotation")
			conn.execute("ALTER TABLE snp_annotation_converted RENAME TO snp_annotation")
			conn.commit()
			if scheme_info:
				write_scheme_info(conn, scheme_info)
			conn.execute("VACUUM")
		finally:
			conn.close()
		try:
			check_database(partial)
		except ConnectionError as e:
			raise ConversionError("The converted copy of {db} failed the database check, nothing was written: {error}".format(db=source_db, error=e))
		os.replace(partial, target_db)
	except sqlite3.Error as e:
		raise ConversionError("The database {db} could not be converted, nothing was written: {error}".format(db=source_db, error=e))
	finally:
		_remove(partial)
	if sequences is not None:
		sequences.update(placed)
	return len(converted_rows)


def _joined_rows(path):
	conn = sqlite3.connect(path)
	try:
		columns = set(row[1] for row in conn.execute("PRAGMA table_info(snp_annotation)").fetchall())
		if "sequence_position" in columns:
			raise ConversionError("The database {path} already stores positions per reference sequence, it is in the current format".format(path=path))
		if "position" not in columns:
			raise ConversionError("The database {path} has no position column in snp_annotation".format(path=path))
		return conn.execute('''SELECT snp_annotation.node_id, snp_annotation.snp_id, snp_annotation.position,
			snp_annotation.ancestral_base, snp_annotation.derived_base, snp_annotation.reference,
			snp_annotation.date, snp_annotation.genome_i, snp_references.genome
			FROM snp_annotation LEFT JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)
			ORDER BY snp_annotation.node_id''').fetchall()
	finally:
		conn.close()


def database_genomes(path):
	'''the genomes that carry SNPs in a database'''
	conn = sqlite3.connect(path)
	try:
		return [row[0] for row in conn.execute('''SELECT DISTINCT snp_references.genome FROM snp_annotation
			JOIN snp_references ON (snp_references.id = snp_annotation.genome_i)''').fetchall()]
	finally:
		conn.close()


def table_header(path):
	'''the header row of a SNP or canSNP source file, as a list of column names'''
	with read_text(path) as fin:
		for line in fin:
			line = line.rstrip("\r\n")
			if line.strip() == "" or line.startswith("#"):
				continue
			return line.split("\t")
	return []


def table_stage(path):
	'''What a SNP or canSNP source file needs: "joined" for a position along the joined
		reference sequences (before 2.3.0), "snp" for one canSNP per node without the rank,
		score and typed columns (2.3.0), or "current"'''
	header = table_header(path)
	if "sequence_position" not in header:
		return "joined"
	return "current" if "typed" in header else "snp"


def check_table_format(path):
	'''Raise ConversionError when a SNP source file is already in the current format, before any
		reference is fetched'''
	if table_stage(path) == "current":
		raise ConversionError("{source} already has the rank, score and typed columns, it is in the current "
			"format, nothing to convert".format(source=path))


def database_stage(path):
	'''What a database needs: "joined" for positions along the joined reference sequences
		(before 2.3.0), "snp" for one canSNP per node without the rank, score and typed
		columns (2.3.0), or "current"'''
	conn = sqlite3.connect(path)
	try:
		columns = set(row[1] for row in conn.execute("PRAGMA table_info(snp_annotation)").fetchall())
	finally:
		conn.close()
	if "sequence_position" not in columns:
		return "joined"
	return "current" if "typed" in columns else "snp"


def check_database_format(path):
	'''Raise ConversionError when a database is already in the current format'''
	if database_stage(path) == "current":
		raise ConversionError("The database {path} already stores every canSNP with its rank, score and whether "
			"the scheme types it, it is in the current format, nothing to convert".format(path=path))


def table_genomes(path):
	'''the genomes named in the genome column of a source file'''
	genomes = []
	header = None
	with read_text(path) as fin:
		for line in fin:
			line = line.rstrip("\r\n")
			if line.strip() == "" or line.startswith("#"):
				continue
			fields = line.split("\t")
			if header is None:
				header = fields
				if "genome" not in header:
					raise ConversionError("{path} has no genome column".format(path=path))
				genome_i = header.index("genome")
				continue
			if len(fields) > genome_i and fields[genome_i] not in genomes:
				genomes.append(fields[genome_i])
	return genomes


def references_from_database(path):
	'''{genome: (genbank_id, assembly_name)} from the snp_references table of a database'''
	conn = sqlite3.connect(path)
	try:
		return {genome: (genbank_id, assembly) for genome, genbank_id, assembly in
			conn.execute("SELECT genome, genbank_id, assembly_name FROM snp_references").fetchall()}
	finally:
		conn.close()


def references_from_file(path):
	'''{genome: (genbank_id, assembly_name)} from a references source file'''
	references = {}
	header = None
	with read_text(path) as fin:
		for line in fin:
			line = line.rstrip("\r\n")
			if line.strip() == "" or line.startswith("#"):
				continue
			fields = line.split("\t")
			if header is None:
				header = fields
				missing = [c for c in ("genome", "genbank_id", "assembly_name") if c not in header]
				if missing:
					raise ConversionError("{path} has no column(s) {missing}".format(path=path, missing=", ".join(missing)))
				continue
			row = dict(zip(header, fields))
			references[row["genome"]] = (row.get("genbank_id"), row.get("assembly_name"))
	return references


class _Results(list):
	'''collects what DownloadGenomes.download reports, without a process queue'''
	def put(self, item):
		self.append(item)


def genome_files(refdir, genomes, references=None):
	'''{genome: path} of the reference files <refdir>/<genome>.fna; a genome without its file is
		downloaded from NCBI GenBank with the accessions of references
		({genome: (genbank_id, assembly_name)}, from snp_references or a references file)'''
	files = {}
	missing = []
	for genome in genomes:
		if genome in files:
			continue
		path = os.path.join(refdir, "{genome}.fna".format(genome=genome))
		if os.path.exists(path):
			files[genome] = path
		else:
			missing.append(genome)
	if not missing:
		return files
	references = references or {}
	unknown = [genome for genome in missing if not (references.get(genome) and references[genome][0]
		and references[genome][0] != "NA" and references[genome][1])]
	if unknown:
		raise ConversionError(("The reference file(s) {files} are not in {refdir}, and no GenBank accession and assembly "
			"name are known for {genomes}; place the files in {refdir}").format(
			files=listing("{0}.fna".format(g) for g in unknown), refdir=refdir, genomes=listing(unknown)))
	from CanSNPer2.DownloadGenomes import DownloadGenomes
	downloader = DownloadGenomes(None, source="genbank", directory=refdir)
	failed = []
	for genome in missing:
		genbank_id, assembly = references[genome]
		logger.info("Downloading the reference {genome} ({acc}) into {refdir}".format(genome=genome, acc=genbank_id, refdir=refdir))
		results = _Results()
		downloader.download(genbank_id, assembly, genome, source="genbank", directory=refdir, results=results, snps=None)
		success, reason = (results[0][1], results[0][2]) if results else (False, "the download did not report")
		path = os.path.join(refdir, "{genome}.fna".format(genome=genome))
		if not success or not os.path.exists(path):
			failed.append("{genome} ({acc}): {reason}".format(genome=genome, acc=genbank_id, reason=reason or "no file was written"))
			continue
		files[genome] = path
	if failed:
		raise ConversionError("The reference(s) could not be downloaded: {failed}".format(failed="; ".join(failed)))
	return files
