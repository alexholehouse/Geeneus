[DEVELOPMENT STAGE - NOT FOR GENERAL USE]

Early stage development for a hyper simple API to the NCBI databases.

Entrez provides eUtils for accessing their databases. However, they are not paricualrly easy to get working quickly, and the Biopython.Entrez module is a fairly thin wrapper. Notably, error support is poor (no error status flags) and just the act of getting the sequence for a specific gene is fairly challenging, especially to a non-programming literate user.

Geeneus aims to remove that difficulty, by abstracting away the complications. There may be an atempt to integrate into Biopython once the API is more complete, but for now it runs as a seperate module with Biopython dependencies. 

