from pathlib import Path
from typing import List
from typing import Tuple
from typing import Union

class ProteinSequence:
    def __init__(self, header: str, sequence: str):
        self.header = header
        self.seq_str = sequence

    def __str__(self):
        return f">{self.header}\n{self.seq_str}"

    def __repr__(self):
        return f">{self.header}\n{self.seq_str}"

    def __len__(self):
        return len(self.seq_str)

    def __getitem__(self, index):
        """This should allow for slicing of the sequence"""
        return ProteinSequence(self.header, self.seq_str[index])

    def __add__(self, other):
        """This should allow for concatenation of the sequence"""
        if isinstance(other, ProteinSequence):
            return ProteinSequence(self.header, self.seq_str + other.seq_str)
        elif isinstance(other, str):
            return ProteinSequence(self.header, self.seq_str + other)
        else:
            raise TypeError(f"Cannot concatenate ProteinSequence with {type(other)}")

    def __radd__(self, other):
        """This should allow for concatenation of the sequence"""
        if isinstance(other, ProteinSequence):
            return ProteinSequence(self.header, other.seq_str + self.seq_str)
        elif isinstance(other, str):
            return ProteinSequence(self.header, other + self.seq_str)
        else:
            raise TypeError(f"Cannot concatenate ProteinSequence with {type(other)}")


def parse_header(header_line: str):
    header = header_line.rstrip()
    assert header.startswith(
        ">"
    ), f"Expected header to start with >, instead saw: {header_line}"
    return header[1:]


def import_a3m(filepath):
    with open(filepath) as f:
        info_line = f.readline().rstrip()
        if not info_line.startswith("#"):
            raise ValueError(
                f"{filepath} should begin with #, instead saw: {info_line}"
            )
        lines = [i.strip() for i in f.readlines()]
        query_header = parse_header(lines[0])
        query_sequence = lines[1].rstrip()
        query = ProteinSequence(query_header, query_sequence)
        sequences = []
        for c in range(2, len(lines), 2):
            name = parse_header(lines[c])
            seq = lines[c + 1].rstrip()
            sequences.append(ProteinSequence(name, seq))
    return info_line, query, sequences


def parse_info_line(info_line: str):
    """ """
    info = info_line.split("\t")
    assert (
        len(info) == 2
    ), f"Expected 2 fields separated by a tab in the info line, instead saw: {info_line}"
    query_seq_lengths = info[0][1:].split(",")
    query_seq_cardinality = info[1].split(",")
    return query_seq_lengths, query_seq_cardinality


class MSAa3m:

    def __init__(
        self, info_line: str, query: ProteinSequence, sequences: List[ProteinSequence]
    ):
        self.info_line = info_line
        self.query = query
        assert "-" not in self.query.seq_str, "Query sequence should not contain gaps"
        assert all(
            [i.isupper() for i in self.query.seq_str]
        ), "Query sequence should not contain insertions"
        self.sequences = sequences
        self.query_seq_lengths, self.query_seq_cardinality = self._parse_info_line(
            self.info_line
        )

    def __str__(self):
        msa_str = f"{self.info_line}\n{self.query}\n"
        for seq in self.sequences:
            msa_str += f"{seq}\n"
        return msa_str

    def __repr__(self):
        msa_str = f"{self.info_line}\n{self.query}\n"
        for seq in self.sequences:
            msa_str += f"{seq}\n"
        return msa_str

    @classmethod
    def from_a3m_file(cls, file_path: Union[str, Path]):
        info_line, query, sequences = import_a3m(file_path)
        return cls(info_line, query, sequences)

    @staticmethod
    def _parse_info_line(info_line: str):
        """
        #157,235	1,1
        #398	1
        """
        info = info_line.split("\t")
        assert (
            len(info) == 2
        ), f"Expected 2 fields separated by a tab in the info line, instead saw: {info_line}"
        query_seq_lengths = info[0][1:].split(",")
        query_seq_cardinality = info[1].split(",")
        return query_seq_lengths, query_seq_cardinality

    def _slice_alignment(self, start, end):
        if len(self.query_seq_lengths) > 1:
            raise NotImplementedError(
                "Cannot slice alignment with multiple query sequences yet"
            )
        sliced_query = ProteinSequence(self.query.header, self.query.seq_str[start:end])
        new_info_line = f"#{len(sliced_query.seq_str)}\t{self.query_seq_cardinality[0]}"
        sliced_sequences = []
        for seq in self.sequences:
            core_index = 0
            sliced_seq = ""
            for char in seq.seq_str:
                if core_index >= end:
                    break
                if char.isupper():
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
                elif char.islower():
                    if core_index > start:
                        sliced_seq += char
                elif char == "-":
                    if core_index >= start:
                        sliced_seq += char
                    core_index += 1
            sliced_sequences.append(ProteinSequence(seq.header, sliced_seq))
        return MSAa3m(new_info_line, sliced_query, sliced_sequences)

    def __getitem__(self, index):
        """This should allow for slicing of the sequence"""
        if isinstance(index, slice):
            return self._slice_alignment(index.start, index.stop)
        elif isinstance(index, int):
            return self._slice_alignment(index, index + 1)
        else:
            raise TypeError("Only integers and slices are valid indices for A3M_MSA")

    def __add__(self, other):
        """This should allow for concatenation of the msas."""
        if isinstance(other, MSAa3m):
            new_info_line = f"#{','.join(self.query_seq_lengths)},{','.join(other.query_seq_lengths)}\t{','.join(self.query_seq_cardinality)},{','.join(other.query_seq_cardinality)}"

            # change other query header to a new index.
            # create a new object to avoid changing the original
            last_self_int = [int(i) for i in self.query.header.split("\t")][-1]
            other_ints = [int(i) for i in other.query.header.split("\t")]
            c = last_self_int + 1
            other_int_map = {}
            new_other_int_list = []
            for i in other_ints:
                other_int_map[str(i)] = str(c)
                new_other_int_list.append(str(c))
                c += 1
            new_other_header = "\t".join(new_other_int_list)
            other_query = ProteinSequence(new_other_header, other.query.seq_str)

            new_query = ProteinSequence(
                f"{self.query.header}\t{other_query.header}",
                f"{self.query.seq_str}{other_query.seq_str}",
            )
            new_sequences = []
            if len(self.query_seq_lengths) == 1:
                new_sequences.append(self.query + "-" * len(other_query))
            for seq in self.sequences:
                new_seq = seq + "-" * len(other_query)
                if not all([i == "-" for i in new_seq.seq_str]):
                    new_sequences.append(new_seq)
            if len(other.query_seq_lengths) == 1:
                new_sequences.append("-" * len(self.query) + other_query)
            for seq in other.sequences:
                new_seq = "-" * len(self.query) + seq
                if new_seq.header in other_int_map.keys():
                    new_seq.header = other_int_map[new_seq.header]
                if not all([i == "-" for i in new_seq.seq_str]):
                    new_sequences.append(new_seq)
            # add the original 2 queries back to the MSA but only
            # if the MSA was not previously concatenated
            return MSAa3m(new_info_line, new_query, new_sequences)
        else:
            raise TypeError(f"Cannot concatenate A3M_MSA with {type(other)}")

    def save(self, filepath):
        with open(filepath, "w") as f:
            f.write(str(self))
