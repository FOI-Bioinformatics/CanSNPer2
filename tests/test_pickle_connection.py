import multiprocessing
import pickle
import sqlite3

from CanSNPer2.DownloadGenomes import DownloadGenomes
from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions


def make_db(path):
	conn = sqlite3.connect(path)
	conn.executescript("""
		CREATE TABLE nodes (id integer PRIMARY KEY, name text NOT NULL);
		CREATE TABLE snp_references (id integer PRIMARY KEY, genome varchar(6), strain varchar(200),
			genbank_id varchar(20), refseq_id varchar(20), assembly_name varchar(200));
		INSERT INTO nodes (id, name) VALUES (1, 'N1'), (2, 'N2');
	""")
	conn.commit()
	conn.close()


def count_nodes(db, queue):
	queue.put(db.query("SELECT count(*) FROM nodes").fetchone()[0])


def test_connection_object_can_be_pickled(tmp_path):
	path = tmp_path / "test.db"
	make_db(path)
	db = CanSNPdbFunctions(str(path))
	copy = pickle.loads(pickle.dumps(db))
	assert copy.conn is not db.conn
	assert copy.query("SELECT name FROM nodes ORDER BY id").fetchall() == [("N1",), ("N2",)]
	db.disconnect()
	copy.disconnect()


def test_connection_object_reaches_a_spawned_worker(tmp_path):
	'''spawn is the default start method on macOS, and forkserver, which also pickles the
		worker arguments, is the default on linux from Python 3.14'''
	path = tmp_path / "test.db"
	make_db(path)
	db = CanSNPdbFunctions(str(path))
	ctx = multiprocessing.get_context("spawn")
	queue = ctx.Queue()
	worker = ctx.Process(target=count_nodes, args=(db, queue))
	worker.start()
	result = queue.get(timeout=60)
	worker.join(timeout=60)
	db.disconnect()
	assert worker.exitcode == 0
	assert result == 2


def test_download_worker_target_can_be_pickled(tmp_path):
	'''DownloadGenomes.run starts a Process with the bound method self.download as target,
		so the whole object, including its database connection, is pickled under spawn'''
	path = tmp_path / "test.db"
	make_db(path)
	downloader = DownloadGenomes(str(path), directory=str(tmp_path / "references"))
	copy = pickle.loads(pickle.dumps(downloader.download))
	assert copy.__self__.database.query("SELECT count(*) FROM nodes").fetchone()[0] == 2
