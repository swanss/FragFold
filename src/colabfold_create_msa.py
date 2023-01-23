def readFastaLines(file):
    header = file.readline().rstrip()
    if header != '' and not header.startswith('>'):
        raise ValueError('Header line in FASTA should begin with >, instead saw:'+header)
    sequence = file.readline().rstrip()
    return header,sequence

def extractSubsequence(sequence,seqRange:tuple):
    # Remember: residue number is 1-indexed, but the sequence is zero-indexed
    return sequence[seqRange[0]-1:seqRange[1]]

def hasGaps(sequence,gapchar:str='-'):
    # search for gaps within the sequence (between non-gap chars)
    return sum([x for x in map(lambda x: 1 if x != '' else 0,sequence.split(gapchar))]) > 1

def hasLower(sequence):
    for a in sequence:
        if a.islower():
            return True
    return False

def countLower(sequence):
    count = 0
    for a in sequence:
        if a.islower():
            count+=1
    return count

def createMSA(inputMSA: str, protein_range: tuple, fragment_range: tuple, subsample: int, a3m_output_path: str, protein_copies: int = 1):
    full_query_seq = ''
    new_info_line = ''

    with open(inputMSA,'r') as input_file, open(a3m_output_path,'w') as output_file:
        info_line = input_file.readline().rstrip()
        # print('info line:',info_line)
        if not info_line.startswith('#'):
            raise ValueError('File should begin with #, instead saw'+info_line)
            
        # First line is whole protein sequence
        header,whole_sequence = readFastaLines(input_file)
        query_protein_sequence = extractSubsequence(whole_sequence,protein_range)
        query_fragment_sequence = extractSubsequence(whole_sequence,fragment_range)
        full_query_seq = query_protein_sequence + query_fragment_sequence
        new_info_line = f"#{len(query_protein_sequence)},{len(query_fragment_sequence)}\t{protein_copies},1"+'\n'
        
        output_file.write(new_info_line)
        output_file.write('>101\t102\n')
        output_file.write(full_query_seq+'\n')
            
        count = 1
        # Remaining lines are the matches from the MSA
        while header != '':
            header,sequence = readFastaLines(input_file)
            if subsample > 0 and count == subsample:
                break
            
            # Get the protein/fragment sequence
            protein_sequence = extractSubsequence(sequence,protein_range)
            fragment_sequence = extractSubsequence(sequence,fragment_range)

            if protein_sequence.replace('-','') != '':
                new_header = header + '\t101\n'
                new_protein_sequence = protein_sequence.ljust(len(full_query_seq),'-') 
                new_protein_sequence = new_protein_sequence + '-'*countLower(protein_sequence) + '\n'
                output_file.write(new_header)
                output_file.write(new_protein_sequence)
            if fragment_sequence.replace('-','') != '' and not hasGaps(fragment_sequence):
                new_header = header + '\t102\n'
                new_fragment_sequence = fragment_sequence.rjust(len(full_query_seq),'-')
                new_fragment_sequence = new_fragment_sequence + '-'*countLower(fragment_sequence) + '\n'
                output_file.write(new_header)
                output_file.write(new_fragment_sequence)
            count+=1

def createMSAHeteromicInteraction(fullLengthProteinMSA: str, protein_range: tuple, fragmentProteinMSA: str, fragment_range: tuple, subsample: int, a3m_output_path: str, protein_copies: int = 1):
    full_query_seq = ''
    new_info_line = ''

    with open(fullLengthProteinMSA,'r') as fulllength_input_file, open(fragmentProteinMSA,'r') as fragment_input_file, open(a3m_output_path,'w') as output_file:
        info_line = fulllength_input_file.readline().rstrip()
        if not info_line.startswith('#'):
            raise ValueError(f"{fullLengthProteinMSA} should begin with #, instead saw'+info_line")
        info_line = fragment_input_file.readline().rstrip()
        if not info_line.startswith('#'):
            raise ValueError(f"{fragmentProteinMSA} should begin with #, instead saw'+info_line")
            
        # First line in the full-length msa is whole protein sequence
        fulllength_header,whole_sequence = readFastaLines(fulllength_input_file)
        query_protein_sequence = extractSubsequence(whole_sequence,protein_range)

        # First line in fragment msa is the whole sequence of the protein that was broken into fragments
        fragment_header,whole_fragmented_sequence = readFastaLines(fragment_input_file)
        query_fragment_sequence = extractSubsequence(whole_fragmented_sequence,fragment_range)
        full_query_seq = query_protein_sequence + query_fragment_sequence
        new_info_line = f"#{len(query_protein_sequence)},{len(query_fragment_sequence)}\t{protein_copies},1"+'\n'
        
        output_file.write(new_info_line)
        output_file.write('>101\t102\n')
        output_file.write(full_query_seq+'\n')
            
        # Remaining lines are the matches from the MSAs
        # First get full-length MSA entries
        count = 1
        while fulllength_header != '':
            fulllength_header,sequence = readFastaLines(fulllength_input_file)
            if subsample > 0 and count == subsample:
                break
            # Get the protein/fragment sequence
            protein_sequence = extractSubsequence(sequence,protein_range)
            if protein_sequence.replace('-','') != '':
                new_header = fulllength_header + '\t101\n'
                new_protein_sequence = protein_sequence.ljust(len(full_query_seq),'-') 
                new_protein_sequence = new_protein_sequence + '-'*countLower(protein_sequence) + '\n'
                output_file.write(new_header)
                output_file.write(new_protein_sequence)
            count+=1
        # Next get fragment MSA entries
        count = 1
        while fragment_header != '':
            fragment_header,sequence = readFastaLines(fragment_input_file)
            if subsample > 0 and count == subsample:
                break
            # Get the protein/fragment sequence
            fragment_sequence = extractSubsequence(sequence,fragment_range)
            if fragment_sequence.replace('-','') != '' and not hasGaps(fragment_sequence):
                new_header = fragment_header + '\t102\n'
                new_fragment_sequence = fragment_sequence.rjust(len(full_query_seq),'-')
                new_fragment_sequence = new_fragment_sequence + '-'*countLower(fragment_sequence) + '\n'
                output_file.write(new_header)
                output_file.write(new_fragment_sequence)
            count+=1
        
        
