import unittest
import Geeneus


class TestProteomeFunctions(unittest.TestCase):

    # Build manager object for all tests here
    def setUp(self):
        self.manager = Geeneus.Proteome.ProteinManager("alex.holehouse@gmail.com")

    def test_get_sequence(self):
        # make sure we can pull a protein sequence

        testSeq = "maaaaaagagpemvrgqvfdvgprytnlsyigegaygmvcsaydnvnkvrvaikkispfehqtycqrtlreikillrfrheniigindiiraptieqmkdvyivqdlmetdlykllktqhlsndhicyflyqilrglkyihsanvlhrdlkpsnlllnttcdlkicdfglarvadpdhdhtgflteyvatrwyrapeimlnskgytksidiwsvgcilaemlsnrpifpgkhyldqlnhilgilgspsqedlnciinlkarnyllslphknkvpwnrlfpnadskaldlldkmltfnphkrieveqalahpyleqyydpsdepiaeapfkfdmelddlpkeklkelifeetarfqpgyrs"
        
        self.assertEqual(testSeq, self.manager.get_protein_sequence("NP_002736"))
        
        
if __name__ == '__main__':
    print "Please run all tests together"
