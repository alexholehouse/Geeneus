import unittest
from geeneus import backend

class TestProteinParserFunctions(unittest.TestCase):

    def setUp(self):
        
        # build parser object with cache turned off
        self.ParserObject = backend.ProteinParser.ProteinRequestParser("alex.holehouse@gmail.com", cache=True, retry=4)

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
        print "Size = " + str(self.ParserObject.get_size_of_datastore())
        IDLIST = ["AAB29246", "1842230", "Q13480","NP_997006","P16144","NP_000204","NP_001005619","gi|61742777","Q9UIQ6","NP_005566","NP_787116","gi|8923710","Q9UKY7","NP_060018","gi|20986529","NP_002736","gi|4503787","P42685","NP_002022", "gi|20562757","B0LPG3","B3KR49","P27361","Q7Z3H5","NP_001035145","NP_002737","XP_055766","IPI00018195","gi|20357552","Q14247","Q53HG7","Q76MU0","NP_005222","IPI00029601","gi|21361602","A8K6H9","Q96QD8","Q96QD8-1","Q9HAV3","NP_061849","IPI00410034","P49023","P49023-1","PAXI_HUMAN","NP_002850","IPI00335634","C9JH84","NP_065939","IPI00001632","IPI00641889","gi|4759082","A8K2W3","O95810","SDPR_HUMAN","NP_004648","IPI00005809","gi|9257199","Q9UQB8","NP_001138360","NP_059344","NP_059345","gi|37551295","gi|50843820","gi|51468542","Q5T5P2","NP_001091970","NP_062536","XP_166112","XP_495798","gi|27886584","Q14289","NP_004094","NP_775266","NP_775268","gi|4505055","A8K379","P07948","P07948-1","NP_001104567","NP_002341","IPI00298625","gi|50845388","P07355-2","NP_001002858","IPI00418169","gi|52426745","P22681","NP_005179","IPI00027269","Q6ZVM7","Q6ZVM7-1","NP_001076437","IPI00446294","P30530","P30530-1","NP_068713","IPI00296992","Q01804","Q01804-1","NP_059963","IPI00399254","Q86YV5","NP_001074295","IPI00166578","IPI00739386","gi|15451856","A9XTE5","Q03135","Q03135-1","Q2TNI1","NP_001744","IPI00009236","gi|5803137","P98179","NP_006734","IPI00024320","gi|4758248","P98172","NP_004420","IPI00024307","gi|21361181","P05023","NP_000692","NP_001153705","gi|16936528","P24941","NP_001789","gi|4506403","A8K556","B3KV45","Q8NFJ5","NP_003970","IPI00022624"]
        

        outlist = self.ParserObject.batchFetch(self.ParserObject.get_sequence, IDLIST)
        print "Size = " + str(self.ParserObject.get_size_of_datastore())


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
        
