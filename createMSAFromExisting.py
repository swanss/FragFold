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

def createMSA(inputMSA: str, protein_range: tuple, fragment_range: tuple, subsample: int, a3m_output_path: str):
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
        new_info_line = f"#{len(query_protein_sequence)},{len(query_fragment_sequence)}\t1,1"+'\n'
        
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
