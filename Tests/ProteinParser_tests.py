import unittest
import Geeneus.Backend

class TestProteinParserFunctions(unittest.TestCase):

    def setUp(self):
        
        # build parser object with cache turned off
        self.ParserObject = Geeneus.Backend.ProteinParser.ProteinRequestParser("alex.holehouse@gmail.com", False)

        self.assertTrue(not self.ParserObject.error()) 

        self.seqAAB29246 = "msepagdvrqnpcgskacrrlfgpvdseqlsrdcdalmagciqearerwnfdfvtetplegdfawervrglglpklylptgprrgrdelgggrrpgtspallqgtaeedhvdlslsctlvprsgeqaegspggpgdsqgrkrrqtsmtdfyhskrrlifskrkp"
        self.seq1842230 = "meepqsdlsielplsqetfsdlwkllppnnvlstlpssdsieelflsenvtgwledsggalqgvaaaaastaedpvtetpapvasapatpwplsssvpsyktfqgdygfrlgflhsgtaksvtctyspslnklfcqlaktcpvqlwvnstpppgtrvramaiykklqymtevvrrcphherssegdslappqhlirvegnlhaeylddkqtfrhsvvvpyeppevgsdcttihynymcnsscmggmnrrpiltiitledpsgnllgrnsfevricacpgrdrrteeknfqkkgepcpelppksakralptntssspppkkktldgeyftlkirgherfkmfqelnealelkdaqaskgsedngahssylkskkgqsasrlkklmikregpdsd"


    def test_get_protein_sequence(self):
                
        self.assertEqual(self.seqAAB29246, self.ParserObject.get_sequence("AAB29246"))

        
    def test_get_fake_protein_sequence(self):
        print"\nYou should see a warning regarding 'FAKETEST' below"
        self.assertEqual("", self.ParserObject.get_sequence("FAKETEST"))


    def test_get_batch_protein_sequence(self):
        print"\nYou should see a warning regarding 'FAKETEST' below"
        IDLIST = ["AAB29246", "FAKETEST", "1842230"]
        outlist = self.ParserObject.batchFetch(self.ParserObject.get_sequence, IDLIST)
            
        self.assertEquals(outlist[IDLIST[0]], self.seqAAB29246)
        self.assertEquals(outlist[IDLIST[1]], "")
        self.assertEquals(outlist[IDLIST[2]], self.seq1842230)


    def test_get_variants(self):

        variants = self.ParserObject.get_variants("P16144")

        if len(variants) > 0:
            self.assertEquals('38', variants[0]['Location'])

        else:
            print "If we're connected to the internet this test has failed!"

    def test_purge(self):
        self.ParserObject.purgeDataStore()
        
