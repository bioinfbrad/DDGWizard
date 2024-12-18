.. _Generate Feature-Enriched ΔΔG Data:

.. raw:: html

    <style>.highlight {
            background-color: #E7FC9F;
            color: #000000;
            padding: 6px;
            font-size: 12px;
            font-weight: 100;
        }
            .keyword-highlight {
            background-color: #FFFFF0;
            color: #FF3366;
            padding: 6px;
            font-size: 12px;
            font-weight: 100;
        }
    </style>

Generate Feature-Enriched ΔΔG Data
====================================

.. raw:: html

    <div style="text-align: justify;">
    This guide is intended to show users how to use DDGWizard to process raw ΔΔG data and output feature-enriched data. It can help user obtain more diverse feature information for their own ΔΔG dataset, facilitating further analysis, feature selection, and machine learning.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>1. Prepare a Blast database</h4>
    <p>A Blast database is required for the program to run. The program will use the path to the Blast database to invoke it and perform sequence alignment.</p>
    <p></p>
    To construct a Blast database, you first need to prepare a <span class="keyword-highlight">fasta</span> file of the protein sequence database.
    <p></p>
    The richness of the sequence database should be abundant. You can use your own <span class="keyword-highlight">fasta</span> database file, but we recommend downloading it from <a href="https://ftp.uniprot.org/pub/databases/uniprot/uniref/">Uniref Databases</a>.
    <p></p>
    Our program was tested using <a href="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/">Uniref50 database</a>.
    <p></p>
    If you download <span class="keyword-highlight">Uniref50 database</span>, you need to unzip it (taking Uniref50 as an example):
    </div>

.. code-block::

    $ gzip -d uniref50.fasta.gz

.. raw:: html

    <div style="text-align: justify;">
    You can obtain a <span class="keyword-highlight">fasta</span> file. Then use Blast suite to create a Blast database using obtained <span class="keyword-highlight">fasta</span> file.
    <p></p>
    Downloaded <span class="keyword-highlight">blast+ 2.13.0</span> should be in the path <span class="keyword-highlight">DDGWizard/bin/ncbi_blast_2_13_0+/</span>. Please use the command as follows:
    <p></p>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cd DDGWizard/bin/ncbi_blast_2_13_0+/bin/
    $ ./makeblastdb -in <b>&lt;the path to fasta file&gt;</b> -dbtype prot -out <b>&lt;the path to save Blast database&gt;</b>/<b>&lt;the name to assign for Blast database&gt;</b> -parse_seqids
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    This step will take some time, depending on the size of the database file and the performance of your computer system.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>2. Running template</h4>
    <p>We first provide you with a running template of running DDGWizard's feature calculation pipeline, and then explain the specifics of each parameter in detail.</p>
    <p></p>
    You can run the program with:
    <p></p>
    <div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ conda activate DDGWizard
    $ cd DDGWizard/
    $ python Generate_Dataset_Executable.py \
        --raw_dataset_path <b>&lt;the path to csv file of raw data&gt;</b> \
        --db_folder_path <b>&lt;the path to save Blast database&gt;</b> \
        --db_name <b>&lt;the name to assign for Blast database&gt;</b> \
        --if_reversed_data 1 \
        --blast_process_num 4 \
        --mode whole \
        --process_num 4 \
        --container_type <b>&lt;Docker or Singularity or None&gt;</b>
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>3. Parameter details</h4>
    Below are the details of the parameters for program to generate complete  ΔΔG feature set:
    <p></p>
    (1). <span class="keyword-highlight">raw_dataset_path</span>
    <p></p>
    This parameter indicates that you need to provide the path to a <span class="keyword-highlight">csv</span> file, which contains the raw data you want to use to generate ΔΔG feature set.
    <p></p>
    In the path <span class="keyword-highlight">DDGWizard/src</span>, there is a sample file <span class="keyword-highlight">Sample.csv</span> that you can use directly for testing and as a reference.
    <p></p>
    We list some of the contents of this file here, and provide detailed descriptions of each column's attributes in the table file:
    <p></p>
    <div>

+-------------+---------------------------+--------------------+----------------+----------+------------------+
| PDB         | Amino Acid Substitution   | Chain ID           | ddG            |   pH     |  T               |
+=============+===========================+====================+================+==========+==================+
| 1AAR        | K6E                       | A                  | 0.53           |   5      |  25              |
+-------------+---------------------------+--------------------+----------------+----------+------------------+
| 1AAR        | K6Q                       | A                  | 0.26           |   5      |  25              |
+-------------+---------------------------+--------------------+----------------+----------+------------------+
| 1AAR        | H68E                      | A                  | 0.77           |   5      |  25              |
+-------------+---------------------------+--------------------+----------------+----------+------------------+
| ...         | ...                       | ...                |   ...          |  ...     |  ...             |
+-------------+---------------------------+--------------------+----------------+----------+------------------+

.. raw:: html

    <div style="text-align: justify;">
    Description of attributes for each column in the table file:
    <div style="margin-left: 40px;">
    <p></p>
    a. <span class="keyword-highlight">PDB</span>: This attribute requires to provide a <span class="keyword-highlight">PDB</span> identifier sourced from <a href="https://www.rcsb.org/">the RCSB database</a>. Using the <span class="keyword-highlight">PDB</span> identifier program can automatically download the <span class="keyword-highlight">PDB</span> file.
    <p></p>
    b. <span class="keyword-highlight">Amino Acid Substitution</span>: It consists of one-letter code of the wild-type amino acid, the sequential number of the mutation site, and the code of the mutant amino acid, for describing substitution of amino acids caused by the mutation. For example, K6Q represents a substitution where lysine at the 6th position of protein sequence is substituted with glutamine.
    <p></p>
    c. <span class="keyword-highlight">Chain ID</span>: Indicate the protein chain where the mutation site is located.
    <p></p>
    d. <span class="keyword-highlight">ddG</span>: Require to provide the ΔΔG values of users' own raw dataset. For users with machine learning needs, this value can serve as the regression target. If users only require generating features, this attribute can be set to any numerical value without affecting the generation of other features.
    <p></p>
    e. <span class="keyword-highlight">pH</span>: Specify at which pH the mutation occurs.
    <p></p>
    f. <span class="keyword-highlight">T</span>: Specify at which temperature the mutation occurs.
    <p></p>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (2). <span class="keyword-highlight">--db_folder_path</span>
    This parameter indicates the folder path of the Blast database that user have prepared.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (3). <span class="keyword-highlight">--db_name</span>
    This parameter indicates the name of the Blast database that user have prepared.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (4). <span class="keyword-highlight">--if_reversed_data</span>
    This parameter requires user to provide a value of 0 or 1. The value of 0 means only generating features for the direct mutation, while the value of 1 means also generating the features for the reverse mutations of the mutations provided.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (5). <span class="keyword-highlight">--blast_process_num</span>
    This parameter requires user to provide an integer greater than 0 and less than 200. It represents the number of processes (multiprocessing) DDGWizard will use for sequence alignment.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (6). <span class="keyword-highlight">--mode</span>
    Please provide the default value <span class="keyword-highlight">whole</span>.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (7). <span class="keyword-highlight">--process_num</span>
    This parameter requires user to provide an integer greater than 0 and less than 200. It represents the number of processes (multiprocessing) DDGWizard will use for generating features.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (8). <span class="keyword-highlight">--container_type</span>
    This parameter requires user to provide a value of <span class="keyword-highlight">D</span> or <span class="keyword-highlight">S</span> or <span class="keyword-highlight">-</span> (default). The value of <span class="keyword-highlight">D</span> means using <span class="keyword-highlight">Docker</span> as container system, the value of <span class="keyword-highlight">S</span> means using <span class="keyword-highlight">Singularity</span> as container system, and the value of <span class="keyword-highlight">-</span> means skipping running PROFbval.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>4. Output</h4>
    There will be an output <span class="keyword-highlight">csv</span> file <span class="keyword-highlight">features_table.csv</span> located in <span class="keyword-highlight">DDGWizard/src/Feature_Res/</span>, which will record complete generated features.
    <p></p>
    </div>


.. raw:: html

    <div style="text-align: justify;">
    <h4>5. Notes</h4>
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    (1). When running DDGWizard, you need to <span class="keyword-highlight">cd</span> to the top-level directory of the program to execute the program.
    <p></p>
    (2). DDGWizard supports multi-process handling itself. If you wish to run multiple instances of DDGWizard to fully utilize your computer's resources, we recommend using the multi-process parameters provided by DDGWizard.
    <p></p>
    We don't recommend to achieve multi-process handling of DDGWizard by user themselves.
    <p></p>
    If user need to run multiple instances of DDGWizard at the same time by themselves, please avoid running multiple instances of DDGWizard from the same folder, as the program synchronizes files within the folder, which can cause synchronization errors. <b>Please make multiple copies of the DDGWizard folder and run each instance separately in its own folder.</b>
    <p></p>
    (3). <b>Do not place your files in the top-level folder of DDGWizard.</b> DDGWizard will automatically clean files in the top-level folder to maintain multi-process synchronization.
    <p></p>
    (4). <b>The complete log file is saved at the path <span class="keyword-highlight">DDGWizard/src/log.txt</span>.</b>
    <p></p>
    </div>


