--
-- TARGETSCAN_DB schema
--
--

-- Metadata tables.
CREATE TABLE metadata (
  name VARCHAR(80) PRIMARY KEY,
  value VARCHAR(255)
);

CREATE TABLE map_counts (
  map_name VARCHAR(80) PRIMARY KEY,
  count INTEGER NOT NULL
);

CREATE TABLE map_metadata (
  map_name VARCHAR(80) NOT NULL,
  source_name VARCHAR(80) NOT NULL,
  source_url VARCHAR(255) NOT NULL,
  source_date VARCHAR(20) NOT NULL
);

CREATE TABLE genes (
  _id INTEGER PRIMARY KEY,
  gene_id VARCHAR(10) NOT NULL UNIQUE                  -- Entrez Gene ID
);

CREATE TABLE mirna_family (
  _id INTEGER PRIMARY KEY,
  name VARCHAR(255) NOT NULL  	                       -- miRNA family
);

CREATE TABLE species (
  id VARCHAR(10) PRIMARY KEY,	                       -- species ID from NCBI
  name VARCHAR(100) NOT NULL                           -- species name
);

CREATE TABLE seed_match (
  _id INTEGER PRIMARY KEY,
  name VARCHAR(10)
);

CREATE TABLE mirbase (
  mirbase_id VARCHAR(50) PRIMARY KEY,                  -- MiRBase ID
  mirbase_accession CHAR(12) NOT NULL UNIQUE,          -- MiRBase accession
  family INTEGER NOT NULL,                             -- REFERENCES family
  seed_m8 CHAR(7) NOT NULL,                            -- seed m8
  species VARCHAR(10) NOT NULL,                        -- REFERENCES species
  mature_sequence VARCHAR(100) NOT NULL,               -- mature sequence
  family_conservation VARCHAR(3) NOT NULL,             -- family convervation
  FOREIGN KEY (family) REFERENCES mirna_family (_id),
  FOREIGN KEY (species) REFERENCES species (id)
);

CREATE TABLE targets (
  family INTEGER NOT NULL,                             -- REFERENCES family
  target INTEGER NOT NULL,                             -- REFERENCES genes
  species VARCHAR(10) NOT NULL,                        -- REFERENCES species
  utr_start INTEGER NOT NULL,                          -- UTR start
  utr_end INTEGER NOT NULL,                            -- UTR end
  msa_start iNTEGER NOT NULL,                          -- MSA start
  msa_end INTEGER NOT NULL,                            -- MSA end
  seed_match INTEGER NOT NULL,                         -- REFERENCES seed_match
  pct VARCHAR(10) NOT NULL,                            -- PCT
  FOREIGN KEY (family) REFERENCES mirna_family (_id),
  FOREIGN KEY (target) REFERENCES genes (_id),
  FOREIGN KEY (species) REFERENCES species (id),
  FOREIGN KEY (seed_match) REFERENCES seed_match (_id)
);

-- Indexes

