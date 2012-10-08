import unittest
import Geeneus.Backend

class TestProteinParserFunctions(unittest.TestCase):
    
    def setUp(self):
        
        # build parser object with cache turned off
        self.ParserObject = Geeneus.Backend.ProteinParser.ProteinRequestParser("alex.holehouse@gmail.com", False)

        self.assertTrue(not self.ParserObject.error()) 

    def test_get_protein_sequence(self):
        self.assertEqual("msepagdvrqnpcgskacrrlfgpvdseqlsrdcdalmagciqearerwnfdfvtetplegdfawervrglglpklylptgprrgrdelgggrrpgtspallqgtaeedhvdlslsctlvprsgeqaegspggpgdsqgrkrrqtsmtdfyhskrrlifskrkp",self.ParserObject.get_sequence("AAB29246"))
        

    def test_get_fake_protein_sequence(self):
        self.assertEqual("", self.ParserObject.get_sequence("FAKETEST"))

    def test_get_batch_protein_sequence(self):
        IDLIST = ["AAB29246", "FAKETEST", "1842230"]
        outlist = self.ParserObject.batchFetch(self.ParserObject.get_sequence, IDLIST)
            
        self.assertEquals(outlist[IDLIST[0]],"msepagdvrqnpcgskacrrlfgpvdseqlsrdcdalmagciqearerwnfdfvtetplegdfawervrglglpklylptgprrgrdelgggrrpgtspallqgtaeedhvdlslsctlvprsgeqaegspggpgdsqgrkrrqtsmtdfyhskrrlifskrkp")
        self.assertEquals(outlist[IDLIST[1]], "")
        self.assertEquals(outlist[IDLIST[2]],'meepqsdlsielplsqetfsdlwkllppnnvlstlpssdsieelflsenvtgwledsggalqgvaaaaastaedpvtetpapvasapatpwplsssvpsyktfqgdygfrlgflhsgtaksvtctyspslnklfcqlaktcpvqlwvnstpppgtrvramaiykklqymtevvrrcphherssegdslappqhlirvegnlhaeylddkqtfrhsvvvpyeppevgsdcttihynymcnsscmggmnrrpiltiitledpsgnllgrnsfevricacpgrdrrteeknfqkkgepcpelppksakralptntssspppkkktldgeyftlkirgherfkmfqelnealelkdaqaskgsedngahssylkskkgqsasrlkklmikregpdsd')
