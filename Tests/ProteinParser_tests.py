import unittest
import Geeneus.Backend

class TestProteinParserFunctions(unittest.TestCase):

    def setUp(self):
        
        # build parser object with cache turned off
        self.ParserObject = Geeneus.Backend.ProteinParser.ProteinRequestParser("alex.holehouse@gmail.com", cache=True)

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

    def test_massive_batch_request(self):
        IDLIST = ["AAB29246", "1842230", "Q13480","NP_997006","P16144","NP_000204","NP_001005619","gi|61742777","Q9UIQ6","NP_005566","NP_787116","gi|8923710","Q9UKY7","NP_060018","gi|20986529","NP_002736","gi|4503787","P42685","NP_002022"]

        outlist = self.ParserObject.batchFetch(self.ParserObject.get_sequence, IDLIST)

    def test_get_variants(self):

        variants = self.ParserObject.get_variants("P16144")

        if len(variants) > 0:
            self.assertEquals('38', variants[0]['Location'])

        else:
            print "If we're connected to the internet this test has failed!"

    def test_get_number_of_items(self):
        print "Size = " + str(self.ParserObject.get_size_of_datastore())

    def test_purge(self):
        self.ParserObject.purge_data_store()
        
