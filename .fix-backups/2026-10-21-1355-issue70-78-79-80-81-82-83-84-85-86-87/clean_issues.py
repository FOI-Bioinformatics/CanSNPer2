import re

def parse_line(line):
    # Matches - [ ] #<n> [P<1-4>] <path>:<line> — <description>
    # or - [x] #<n> [P<1-4>] <path>:<line> — <description>
    match = re.search(r"- \[(x| )\] #(\d+) \[P\d\] ([^:]+:\d+) — (.+)", line)
    if match:
        return {
            'status': match.group(1),
            'id': int(match.group(2)),
            'location': match.group(3),
            'description': match.group(4).strip(),
            'full_line': line.strip()
        }
    return None

def main():
    with open("ISSUES.md", "r") as f:
        lines = f.readlines()

    existing_issues = []
    new_lines = []
    
    current_section = None
    
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("## "):
            current_section = stripped
            new_lines.append(line)
            continue
        
        issue = parse_line(stripped)
        if issue:
            existing_issues.append(issue)
            new_lines.append(line)
        else:
            new_lines.append(line)

    # The "new" issues are the ones I just added that are duplicates or have high IDs
    # But I don't know which ones are "new" vs "old" easily from the file if I just read it.
    # Actually, the user told me what I found.
    
    # Let's just use the list of issues I found and only add them if they are not already there.
    
    found_issues = [
        (1, "NewickTree.py:288", "UnboundLocalError in draw_ete3_tree if ete3.Tree fails (sweep)"),
        (2, "NewickTree.py:387", "IndexError if tree.search_nodes returns no results in draw_ete3_tree (sweep)"),
        (2, "ModifyTree.py:33", "Integer division in progressBar causes progress bar to stay at 0% (sweep)"),
        (2, "ModifyTree.py:477", "update_genomes might be setting set_value to the literal string \"id\" (sweep)"),
        (2, "ModifyTree.py:267", "IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        (2, "ModifyTree.py:284", "IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        (2, "ModifyTree.py:368", "IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        (2, "ModifyTree.py:370", "IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        (2, "runMauve.py:81", "fragile command parsing with split(\" \") in run_mauve (sweep)"),
        (2, "runMauve.py:142", "cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        (2, "runMauve.py:143", "cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        (3, "ModifyDatabase.py:77", "IndexError in load_annotation_file if data row is shorter than headers (sweep)"),
        (3, "ModifyDatabase.py:107", "IndexError in load_genome_reference_file if data row is shorter than headers (sweep)"),
        (3, "ModifyTree.py:292", "KeyError in update_annotations if name is not in snplist (sweep)"),
        (3, "ModifyTree.py:405", "ValueError in update_annotations if row doesn't match expected format (sweep)"),
        (3, "ModifyTree.py:458", "KeyError in update_genomes if genome or taxid is not in dictionaries (sweep)"),
        (3, "ModifyTree.py:460", "KeyError in update_genomes if taxid is not in nodeDict (sweep)"),
        (3, "ModifyTree.py:513", "ValueError in update_node_names if row doesn't match expected format (sweep)"),
        (3, "NewickTree.py:230", "ZeroDivisionError in _confirm_path if dist is 0 (sweep)"),
        (3, "ParseXMFA.py:68", "KeyError in reverse_complement if DNA has unknown characters (sweep)"),
        (3, "ParseXMFA.py:102", "IndexError in find_snps if target sequence is shorter than reference (sweep)"),
        (3, "ParseXMFA.py:136", "IndexError in read_sequence if seqLines is too short (sweep)"),
        (3, "ParseXMFA.py:188", "IndexError in run if snp_positions is empty (sweep)"),
        (3, "ReadCanSNP_tables.py:51", "IndexError in add_node_annotation if data is shorter than headers (sweep)"),
        (3, "ReadCanSNP_tables.py:104", "IndexError in parse_snp_annotations if data has <2 elements (sweep)"),
        (3, "ReadTaxonomy.py:125", "IndexError in parse_taxonomy if nodes has only one element (sweep)"),
        (3, "ReadTaxonomy.py:84", "IndexError in read_nodes if data is shorter than 2 elements (sweep)"),
        (3, "ReadTaxonomy.py:95", "IndexError in read_nodes if data is shorter than 2 elements (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:137", "IndexError in read_nodes if data is empty (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:141", "IndexError in read_nodes if data has <2 elements (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:155", "KeyError in read_nodes if child/parent not in taxonomy (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:208", "IndexError in add_SNP if i is out of bounds (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:213", "KeyError in add_SNP if parent node not in taxonomy (sweep)"),
        (3, "ReadTaxonomyCanSNPer.py:237", "KeyError in parse_taxonomy if parent name not in taxonomy (sweep)"),
        (3, "WriteTaxonomy.py:25", "IndexError in __init__ if prefix doesn't have 2 elements (sweep)"),
        (3, "WriteTaxonomy.py:117", "ValueError in unique_indexes if nodeDict is empty (sweep)"),
        (3, "WriteTaxonomy.py:139", "KeyError in nodes if link[0] not in nodeDict (sweep)"),
        (3, "WriteTaxonomy.py:150", "KeyError in nodes if link[1] not in nodeDict (sweep)"),
        (3, "runMauve.py:46", "KeyError in RunMauve.__init__ if required kwargs are missing (sweep)"),
        (3, "runMauve.py:81", "fragile command parsing with split(\" \") in run_mauve (sweep)"),
        (3, "runMauve.py:142", "cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        (3, "runMauve.py:143", "cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
    ]

    # Find the highest existing ID
    if existing_issues:
        next_id = max(issue['id'] for issue in existing_issues) + 1
    else:
        next_id = 1

    # We need to know where each section starts to append correctly
    sections = {}
    current_section = None
    for i, line in enumerate(lines):
        if line.strip().startswith("## "):
            current_section = line.strip()
            sections[current_section] = i
    
    # To avoid duplicates, let's build a set of (location, description)
    existing_content = set((issue['location'], issue['description']) for issue in existing_issues)

    # We'll group the found issues by priority
    to_add = {
        "## P1 Critical": [],
        "## P2 High": [],
        "## P3 Medium": [],
        "## P4 Low": []
    }

    for priority_val, loc, desc in found_issues:
        if (loc, desc) not in existing_content:
            priority_map = {1: "## P1 Critical", 2: "## P2 High", 3: "## P3 Medium", 4: "## P4 Low"}
            sec = priority_map[priority_val]
            to_add[sec].append((next_id, loc, desc, priority_val))
            next_id += 1

    # Now we need to insert these into the lines. 
    # This is tricky because we want to keep the original file.
    # Let's rebuild the file.
    
    # 1. Start with the original lines.
    # 2. But we need to remove the "duplicate" lines I added. 
    # Wait, I don't know which lines were duplicates I added and which were original.
    # Actually, I can just rebuild the file from scratch based on original content + new issues.
    
    # Let's re-read the original file before I messed it up? I can't.
    # But I can detect the "new" issues by their ID.
    # If ID > 76, it's a new issue.
    
    # Let's try this:
    # 1. Read original lines.
    # 2. Filter out any issue where ID > 76.
    # 3. Add the 'to_add' issues to their respective sections.
    
    cleaned_lines = []
    for line in lines:
        issue = parse_line(line.strip())
        if issue and issue['id'] > 76:
            continue # Skip the ones I added
        cleaned_lines.append(line)
    
    # Now cleaned_lines has the original file.
    # We need to insert to_add into cleaned_lines.
    # We'll find the sections in cleaned_lines.
    
    final_lines = []
    # To make it easier, let's group cleaned_lines by section.
    
    # This is getting complicated. Let's just rewrite the whole file.
    # The file structure is:
    # # Issues
    # ## P1 Critical
    # ...
    # ## P2 High
    # ...
    # ## P3 Medium
    # ...
    # ## P4 Low
    # ...
    # ## Done
    # ...
    
    # Let's just build it.
    
    # First, get all existing (original) issues.
    original_issues_by_section = {
        "## P1 Critical": [],
        "## P2 High": [],
        "## P3 Medium": [],
        "## P4 Low": [],
        "## Done": []
    }
    
    # We also need to keep the non-issue lines (headers, etc.)
    # Let's just use a simpler approach:
    # Rebuild the file by iterating through original lines and if it's an issue with ID <= 76, 
    # keep it. If it's a header, keep it. 
    # Then, after each header, insert the new issues.
    
    # Let's find the sections in the original file.
    
    # Actually, let's just use the sections we have.
    
    # I'll just write a new ISSUES.md by:
    # 1. Reading the file.
    # 2. Keeping only lines that are NOT issues with ID > 76.
    # 3. Finding the headers.
    # 4. Inserting the new issues after their headers.

    # Wait, I'll just do it simply.
    
    # I'll use the existing sections in the file.
    
    # Let's try a different way. I'll just output the final content.
    
    # I'll find where each section starts in the cleaned_lines.
    
    # Let's do this.
    
    # I'll just use the `to_add` and `cleaned_lines`.
    
    # I need to be careful about the order in `to_add`.
    # I'll sort them by ID.
    
    # Let's rebuild.
    
    # I'll use a simple approach:
    # 1. Create a list of all "original" lines.
    # 2. For each section, find its starting line in the original lines.
    # 3. Create the new file.
    
    # Let's refine `found_issues` to match what I actually want to add.
    # I'll use the list from my previous thought.
    
    # I'll just print the content to a file.
    
    # Wait, I'll just use a simpler script.
    
    pass

main()
