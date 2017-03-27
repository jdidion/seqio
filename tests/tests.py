from unittest import TestCase, skipIf
from . import *
from seqio import reader
from xphyle.paths import TempDir, PathSpec, FileSpec, PathVar

class Tests(TestCase):
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
    
    def test_reader(self):
        # The following are equivalent
        r1 = reader(
            os.path.join(self.data_dir, 'test*.1.fq.gz'),
            os.path.join(self.data_dir, 'test*.2.fq.gz'))
        self.assertEquals(FileType.FASTQ, r1.file_type)
        self.assertTrue(r1.paired)
        r2 = reader(
            [os.path.join(self.data_dir, f) for f in (
                'testA.1.fq.gz', 'testB.1.fq.gz')],
            [os.path.join(self.data_dir, f) for f in (
                'testA.2.fq.gz', 'testB.2.fq.gz')])
        self.assertEquals(FileType.FASTQ, r2.file_type)
        self.assertTrue(r2.paired)
        r3 = reader(PathSpec(
            self.data_dir,
            FileSpec(
                PathVar('library'),
                PathVar('pair'),
                template='test{library}.{pair}.fq.gz')))
        self.assertEquals(FileType.FASTQ, r3.file_type)
        self.assertTrue(r3.paired)
        
        records1 = list(r1)
        records2 = list(r2)
        records3 = list(r3)
        self.assertListEqual(records1, records2)
        self.assertListEqual(records1, records3)
        
        names = [rec.name for rec in records1]
        self.assertListEqual(['rec1', 'rec2', 'rec3', 'rec4'], names)
        seqs = [rec.seq for rec in records1]
        self.assertListEqual(
            [b'CCTGTGGG', b'AAGACTTG', b'ATCGGTAG', b'CGCCTGCC'], seqs)
        
        names = [rec.name for rec in records2]
        self.assertListEqual(['rec1', 'rec2', 'rec3', 'rec4'], names)
        seqs = [rec.seq for rec in records2]
        self.assertListEqual(
            [b'CTGTAAGT', b'GCGCAGGG', b'AGATCTCG', b'TGCAAGAA'], seqs)
