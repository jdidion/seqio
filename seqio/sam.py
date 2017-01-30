class Sam(SequenceFormat):
    """SAM/BAM/CRAM format files. Paired-end files must be name-sorted. Does not
    support secondary/supplementary reads.
    """
    name = 'sam'
    aliases = ('bam', 'cram')
    delivers_qualities = True
    lib = OptionalDependency('pysam')
    
    def open(path, mode, **kwargs):
        return self.lib.AlignmentFile(path, mode, **kwargs)
    
    def read_record(self, fileobj):
        for record in fileobj:
            if record.is_secondary or record.is_supplementary:
                continue
            return self._create_record(record)
    
    def read_pair(self, fileobj):
        read1 = read2 = None
        for record in fileobj:
            if record.is_secondary or record.is_supplementary:
                continue
            elif record.is_read1 and not read1:
                read1 = self._create_record(record)
            elif record.is_read2 and not read2:
                read2 = self._create_record(record)
            else:
                raise FormatError(
                    "Expected exactly one read1 and one read2 for read pair "
                    "{}".format((read1 or read2).name))
            if read1 and read2:
                break
        if read1.name != read2.name:
            raise FormatError(
                "Consecutive reads {}, {} in paired-end SAM/BAM file do "
                "not have the same name; make sure your file is name-sorted.",
                reads[0].name, reads[1].name)
        return (read1, read2)
    
    def _create_record(self, record):
        return self.sequence_class(
            name=record.query_name,
            sequence=record.query_sequence,
            qualities=''.join(chr(33 + q) for q in record.query_qualities))
    
    def format_record(self, record):
        record = self.lib.AlignedSegment()
        record.query_name = record.name
        record.query_sequence = record.get_full_sequence_str()
        record.flag = 4 # unmapped
        record.query_qualities = pysam.qualitystring_to_array(
            record.get_quality_str())
        return record
    
    def format_pair(self, read1, read2):
        record1 = self.format_record(read1)
        record1.flag = 77 # paired, unmapped, mate unmapped, first
        record2 = self.format_record(read2)
        record2.flag = 141 # paired, unmapped, mate unmapped, second
        return (record1, record2)
