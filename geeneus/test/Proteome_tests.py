import unittest
from geeneus import Proteome

class TestProteomeFunctions(unittest.TestCase):

    # Build manager object for all tests here
    def setUp(self):
        self.manager_cacheOn = Proteome.ProteinManager("alex.holehouse@gmail.com", cache=True)
        self.manager_cacheOff = Proteome.ProteinManager("alex.holehouse@gmail.com", cache=False)


    def test_NCBI_code_coverage(self):
        acc = "729143" # GI for p21 MAPK

        self.manager_cacheOn.get_ID_type(acc)
        self.manager_cacheOn.get_domains(acc)
        self.manager_cacheOn.get_geneID(acc)
        self.manager_cacheOn.get_gene_name(acc)
        self.manager_cacheOn.get_host(acc)
        self.manager_cacheOn.get_isoforms(acc)
        self.manager_cacheOn.get_other_accessions(acc)
        self.manager_cacheOn.get_protein_name(acc)
        self.manager_cacheOn.get_protein_sequence(acc)
        self.manager_cacheOn.get_protein_sequence_length(acc)
        self.manager_cacheOn.get_raw_xml(acc)
        self.manager_cacheOn.get_record_creation_date(acc)
        self.manager_cacheOn.get_record_source(acc)
        self.manager_cacheOn.get_record_version(acc)
        self.manager_cacheOn.get_species(acc)
        self.manager_cacheOn.get_taxonomy(acc)
        self.manager_cacheOn.get_variants(acc)
        
        self.manager_cacheOff.get_ID_type(acc)
        self.manager_cacheOff.get_domains(acc)
        self.manager_cacheOff.get_geneID(acc)
        self.manager_cacheOff.get_gene_name(acc)
        self.manager_cacheOff.get_host(acc)
        self.manager_cacheOff.get_isoforms(acc)
        self.manager_cacheOff.get_other_accessions(acc)
        self.manager_cacheOff.get_protein_name(acc)
        self.manager_cacheOff.get_protein_sequence(acc)
        self.manager_cacheOff.get_protein_sequence_length(acc)
        self.manager_cacheOff.get_raw_xml(acc)
        self.manager_cacheOff.get_record_creation_date(acc)
        self.manager_cacheOff.get_record_source(acc)
        self.manager_cacheOff.get_record_version(acc)
        self.manager_cacheOff.get_species(acc)
        self.manager_cacheOff.get_taxonomy(acc)
        self.manager_cacheOff.get_variants(acc)


    def test_UNIPROT_code_coverage(self):
        acc = "P38936" # GI for p21 MAPK
        
        self.manager_cacheOn.get_ID_type(acc)
        self.manager_cacheOn.get_domains(acc)
        self.manager_cacheOn.get_geneID(acc)
        self.manager_cacheOn.get_gene_name(acc)
        self.manager_cacheOn.get_host(acc)
        self.manager_cacheOn.get_isoforms(acc)
        self.manager_cacheOn.get_other_accessions(acc)
        self.manager_cacheOn.get_protein_name(acc)
        self.manager_cacheOn.get_protein_sequence(acc)
        self.manager_cacheOn.get_protein_sequence_length(acc)
        self.manager_cacheOn.get_raw_xml(acc)
        self.manager_cacheOn.get_record_creation_date(acc)
        self.manager_cacheOn.get_record_source(acc)
        self.manager_cacheOn.get_record_version(acc)
        self.manager_cacheOn.get_species(acc)
        self.manager_cacheOn.get_taxonomy(acc)
        self.manager_cacheOn.get_variants(acc)
        
        self.manager_cacheOff.get_ID_type(acc)
        self.manager_cacheOff.get_domains(acc)
        self.manager_cacheOff.get_geneID(acc)
        self.manager_cacheOff.get_gene_name(acc)
        self.manager_cacheOff.get_host(acc)
        self.manager_cacheOff.get_isoforms(acc)
        self.manager_cacheOff.get_other_accessions(acc)
        self.manager_cacheOff.get_protein_name(acc)
        self.manager_cacheOff.get_protein_sequence(acc)
        self.manager_cacheOff.get_protein_sequence_length(acc)
        self.manager_cacheOff.get_raw_xml(acc)
        self.manager_cacheOff.get_record_creation_date(acc)
        self.manager_cacheOff.get_record_source(acc)
        self.manager_cacheOff.get_record_version(acc)
        self.manager_cacheOff.get_species(acc)
        self.manager_cacheOff.get_taxonomy(acc)
        self.manager_cacheOff.get_variants(acc)



    def test_get_sequence(self):
        # make sure we can pull a protein sequence
     
        testSeq = "maaaaaagagpemvrgqvfdvgprytnlsyigegaygmvcsaydnvnkvrvaikkispfehqtycqrtlreikillrfrheniigindiiraptieqmkdvyivqdlmetdlykllktqhlsndhicyflyqilrglkyihsanvlhrdlkpsnlllnttcdlkicdfglarvadpdhdhtgflteyvatrwyrapeimlnskgytksidiwsvgcilaemlsnrpifpgkhyldqlnhilgilgspsqedlnciinlkarnyllslphknkvpwnrlfpnadskaldlldkmltfnphkrieveqalahpyleqyydpsdepiaeapfkfdmelddlpkeklkelifeetarfqpgyrs"
        
        self.assertEqual(testSeq, self.manager_cacheOff.get_protein_sequence("NP_002736"))
        self.assertEqual(testSeq, self.manager_cacheOn.get_protein_sequence("NP_002736"))
        self.assertEqual(None, self.manager_cacheOff.get_protein_sequence("SHOULDN'T WORK"))
        self.assertEqual(None, self.manager_cacheOn.get_protein_sequence("SHOULDN'T WORK"))

        
    # check the pdb translation is working OK (uses eSearch)
    def test_PDB_translation(self):
        testSeq = "mahhhhhhmaksglrqdpqstaaatvlkraveldsesrypqalvcyqegidlllqvlkgtkdntkrcnlrekiskymdraenikkyldqekedgkyhkqikieenatgfsyeslfreylnetvtevwiedpyirhthqlynflrfcemlikrpckvktihlltsldegieqvqqsrglqeieeslrshgvllevqysssihdreirfnngwmikigrgldyfkkpqsrfslgycdfdlrpchettvdifhkkhtkni"
       
        self.assertEqual(testSeq, self.manager_cacheOff.get_protein_sequence("2YMB_A"))
        self.assertEqual(testSeq, self.manager_cacheOn.get_protein_sequence("2YMB_A"))
       
        print "The following PDB values are generic, so will return a list of GIs which we can't resolve"
        self.assertEqual(None, self.manager_cacheOff.get_protein_sequence("2YMB"))
        self.assertEqual(None, self.manager_cacheOn.get_protein_sequence("2YMB"))

    def test_get_protein_name(self):
        self.assertEqual('p21 [Homo sapiens]', self.manager_cacheOff.get_protein_name("AAB29246"))
        self.assertEqual('p21 [Homo sapiens]', self.manager_cacheOn.get_protein_name("AAB29246"))

    def test_get_raw_xml(self):
        self.assertEqual("AA", self.manager_cacheOff.get_raw_xml("AAB29246")["GBSeq_moltype"])
        self.assertEqual("AA", self.manager_cacheOn.get_raw_xml("AAB29246")["GBSeq_moltype"])

    def test_get_variants(self):
        self.assertEqual(100, self.manager_cacheOff.get_variants("P42685")[0]["location"])
        self.assertEqual(100, self.manager_cacheOn.get_variants("P42685")[0]["location"])
                
    def test_get_gene_name(self):
        self.assertEqual("vpu", self.manager_cacheOff.get_gene_name("P05923"))
        self.assertEqual("vpu", self.manager_cacheOn.get_gene_name("P05923"))

    def test_get_geneID(self):
        self.assertEqual("2444", self.manager_cacheOff.get_geneID("P42685"))
        self.assertEqual("2444", self.manager_cacheOn.get_geneID("P42685"))

    def test_get_host(self):
        self.assertEqual("Homo sapiens (Human)", self.manager_cacheOff.get_host("P05923"))
        self.assertEqual("Homo sapiens (Human)", self.manager_cacheOn.get_host("P05923"))


    
    def test_get_protein_sequence_length(self):
        self.assertEqual(505, self.manager_cacheOff.get_protein_sequence_length("P42685"))
        self.assertEqual(505, self.manager_cacheOn.get_protein_sequence_length("P42685"))
        
    def test_get_ID_type(self):
        self.assertEqual("UniProtKB/Swiss-Prot", self.manager_cacheOff.get_ID_type("P42685")[1])
        self.assertEqual("UniProtKB/Swiss-Prot", self.manager_cacheOn.get_ID_type("P42685")[1])

    def test_run_translation(self):
        self.assertEqual("1169745", self.manager_cacheOff.run_translation("P42685"))
        self.assertEqual("1169745", self.manager_cacheOn.run_translation("P42685"))

    def test_purge(self):
        self.manager_cacheOff.purge()
        self.manager_cacheOn.purge()
        self.assertEqual(0, self.manager_cacheOff.datastore_size())
        self.assertEqual(0, self.manager_cacheOn.datastore_size())

    def test_batch_get_variants(self):
        IDLIST = ["Q99439", "Q86UX7"]
        self.assertEqual(0, len(self.manager_cacheOff.batch_get_variants(IDLIST)["Q99439"]))
        self.assertEqual(0, len(self.manager_cacheOff.batch_get_variants(IDLIST)["Q86UX7"]))
        self.assertEqual(0, len(self.manager_cacheOn.batch_get_variants(IDLIST)["Q99439"]))
        self.assertEqual(0, len(self.manager_cacheOn.batch_get_variants(IDLIST)["Q86UX7"]))
    
    def test_batch_get_protein_name(self):
        IDLIST = ["Q08AM6", "P35582"]
        self.assertEqual(59, len(self.manager_cacheOff.batch_get_protein_name(IDLIST)["Q08AM6"]))
        self.assertEqual(105, len(self.manager_cacheOff.batch_get_protein_name(IDLIST)["P35582"]))
        self.assertEqual(59, len(self.manager_cacheOn.batch_get_protein_name(IDLIST)["Q08AM6"]))
        self.assertEqual(105, len(self.manager_cacheOn.batch_get_protein_name(IDLIST)["P35582"]))

    def test_Uniprot_fallback(self):
        IDLIST = ["P42858", "A6NC57"]
        self.assertEqual("Huntingtin; AltName: Full=Huntington disease protein; Short=HD protein", self.manager_cacheOn.batch_get_protein_name(IDLIST)["P42858"])
        self.assertEqual("Ankyrin repeat domain-containing protein 62", self.manager_cacheOn.batch_get_protein_name(IDLIST)["A6NC57"])
        

    def test_batch_get_protein_sequence(self):
        IDLIST = ["Q6P5D3", "Q14980"]
        self.assertEqual(1388, len(self.manager_cacheOff.batch_get_protein_sequence(IDLIST)["Q6P5D3"]))
        self.assertEqual(2115, len(self.manager_cacheOff.batch_get_protein_sequence(IDLIST)["Q14980"]))
        self.assertEqual(1388, len(self.manager_cacheOn.batch_get_protein_sequence(IDLIST)["Q6P5D3"]))
        self.assertEqual(2115, len(self.manager_cacheOn.batch_get_protein_sequence(IDLIST)["Q14980"]))


    def test_batch_get_edge_cases(self):
        IDLIST = []
        self.assertEqual({}, self.manager_cacheOn.batch_get_protein_sequence(IDLIST))
       
if __name__ == '__main__':
    print "Please run all tests together"
