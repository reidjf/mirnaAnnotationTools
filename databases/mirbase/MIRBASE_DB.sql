--
-- MIRBASE_DB schema
--

-- metadata tables
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

-- data tables
CREATE TABLE literature_references (
  auto_lit INTEGER PRIMARY KEY,             -- reference key
  medline INTEGER NULL,                     -- pubmed id
  title TEXT,                               -- title
  author TEXT,                              -- authors
  journal TEXT,                             -- citation
  FOREIGN KEY (auto_lit) REFERENCES mirna_literature_references (auto_lit)
);

CREATE TABLE mirna_2_prefam (
  _id INTEGER NOT NULL PRIMARY KEY ,        -- miRNA key
  auto_prefam INTEGER NOT NULL,             -- miRNA family key
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_chromosome_build (
  _id INTEGER NOT NULL,                     -- miRNA key
  xsome VARCHAR(20) NULL,                   -- chromosome
  contig_start INTEGER NULL,                -- start
  contig_end INTEGER NULL,                  -- end
  strand CHAR(2),                           -- strand
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_cluster (
  mirna_id VARCHAR(40) NOT NULL,            -- miRNA key
  member VARCHAR(40) NOT NULL,              -- miRNA IDs of members of cluster
  cluster INTEGER NOT NULL,                 -- miRNA cluster
  FOREIGN KEY (mirna_id) REFERENCES mirna (mirna_id)
);

CREATE TABLE mirna_context (
  _id INTEGER NOT NULL,                     -- miRNA key
  transcript_id VARCHAR(50) NULL,           -- ensembl transcript ID
  overlap_sense CHAR(2) NULL,               -- overlap strand
  overlap_type VARCHAR(20) NULL,            -- overlap type
  number INTEGER NULL,                      -- exon number (?) FIXME(jfr) check
  transcript_source VARCHAR(50) NULL,       -- db name
  transcript_name VARCHAR(50) NULL,         -- transcript ID
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_database_links (
  _id INTEGER NOT NULL,                     -- miRNA key
  db_id TEXT NOT NULL,                      -- database Name
--  comment TEXT,                             -- ... all NA ... (v.14.0)
  db_link TEXT NOT NULL,                    -- database Accession
  db_secondary TEXT,                        -- secondary database
--  other_params TEXT,                        -- ... all NA ... (v.14.0)
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_hairpin (
  mirna_id VARCHAR(40) NOT NULL,            -- miRNA ID
  hairpin TEXT NOT NULL,                    -- folded stem-loop sequence
  mfe REAL NOT NULL,                        -- mimimun fold energy
  FOREIGN KEY (mirna_id) REFERENCES mirna (mirna_id)
);

CREATE TABLE mirna_literature_references (
  _id INTEGER NOT NULL,                     -- miRNA key
  auto_lit INTEGER NOT NULL,                -- reference key
--  comment TEXT,                             -- ... all NA ... (v.14.0)
  order_added INTEGER NULL,                 -- references order ([1], [2], ...)
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_mature (
  auto_mature INTEGER NOT NULL PRIMARY KEY, -- miRNA mature key
  mature_name VARCHAR(40) NOT NULL,         -- mature ID
  mature_acc VARCHAR(20) NOT NULL,          -- mature Accession
  mature_from VARCHAR(4) NULL,              -- mature Sequence start
  mature_to VARCHAR(4) NULL,                -- mature Sequence end
  evidence TEXT,                            -- experimental evidence
  experiment TEXT,                          -- experiment
  similarity TEXT,                          -- similarity to other miRNAs
  FOREIGN KEY (auto_mature) REFERENCES mirna_pre_mature (auto_mature)
);

CREATE TABLE mirna_prefam (
  auto_prefam INTEGER NOT NULL PRIMARY KEY, -- miRNA family key
  prefam_acc VARCHAR(15) UNIQUE NOT NULL,   -- family Accessiom
  prefam_id VARCHAR(40) UNIQUE NOT NULL,    -- family Name
--  description TEXT,                         -- ... all NA ... (v.14.0)
  FOREIGN KEY (auto_prefam) REFERENCES mirna_2_prefam (auto_prefam)
);

CREATE TABLE mirna_pre_mature (
  _id INTEGER NOT NULL,                     -- miRNA key
  auto_mature INTEGER NOT NULL,             -- miRNA mature key
  FOREIGN KEY (_id) REFERENCES mirna (_id)
);

CREATE TABLE mirna_species (
  auto_species INTEGER NOT NULL,            -- species key
  organism VARCHAR(10) NULL,                -- organism
  division VARCHAR(10) NULL,                -- division
  name VARCHAR(100) NULL,                   -- extended name
  taxonomy VARCHAR(200) NULL,               -- taxonomy
  genome_assembly VARCHAR(15) NULL,         -- genome assembly version
  ensembl_db VARCHAR(50) NULL,              -- ensembl version
  FOREIGN KEY (auto_species) REFERENCES mirna (auto_species)
);

CREATE TABLE mirna (
  _id INTEGER PRIMARY KEY,                  -- miRNA key
  mirna_acc VARCHAR(9) UNIQUE NOT NULL,     -- Accession
  mirna_id VARCHAR(40) NOT NULL,            -- ID
  description VARCHAR(100) NOT NULL,        -- description
  sequence TEXT,                            -- stem-loop sequence
  comment TEXT,                             -- comment
--  auto_species INTEGER NOT NULL             -- species key
  organism VARCHAR(10) NULL                 -- organism
);

-- indexes
