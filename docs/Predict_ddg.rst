.. _Predict ΔΔG:

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

Predict ΔΔG
===========

.. raw:: html

    <div style="text-align: justify;">
    This guide is intended to show users how to use DDGWizard to predict ΔΔG.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>1. Running template</h4>
    <p>We first provide you with a running template of running DDGWizard to predict ΔΔG, and then explain the specifics of each parameter in detail.</p>
    <p></p>
    You can run the program with (predicting ΔΔG also requires the prepared Blast database):
    <p></p>
    <div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ conda activate DDGWizard
    $ cd <b>&lt;/path/to/DDGWizard/&gt;</b>
    $ python Predict_ddG_Executable.py \
        --pred_dataset_path src/Sample_Pred.csv \
        --db_folder_path <b>&lt;/folder/to/save/Blast_database/&gt;</b> \
        --db_name <b>&lt;the_name_to_assign_for_Blast database&gt;</b> \
        --if_reversed_data 0 \
        --blast_process_num 4 \
        --mode whole \
        --process_num 4
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>2. Parameter details</h4>
    Below are the details of the parameters for the ΔΔG prediction program:
    <p></p>
    (1). <span class="keyword-highlight">--pred_dataset_path</span>
    <p></p>
    This parameter indicates that you need to provide the path to a <span class="keyword-highlight">csv</span> file, which contains the mutations' basic information for predicting ΔΔG.
    <p></p>
    In the path <span class="keyword-highlight">DDGWizard/src</span>, there is a sample file <span class="keyword-highlight">Sample_Pred.csv</span> that you can use directly for testing and as a reference.
    <p></p>
    We list some of the contents of this file here, and provide detailed descriptions of each column's attributes in the table file:
    <p></p>
    <div>

+-------------+----------------------------+-------------------+----------+----------+
| PDB         | Amino Acid Substitution    | Chain ID          |   pH     |  T       |
+=============+============================+===================+==========+==========+
| 1SHG        | Y57H                       |   A               |   7      |  25      |
+-------------+----------------------------+-------------------+----------+----------+
| 2AFG        | C117I                      |   A               |   7      |  25      |
+-------------+----------------------------+-------------------+----------+----------+
| 2LZM        | M102L                      |   A               |   3      |  52      |
+-------------+----------------------------+-------------------+----------+----------+
| ...         | ...                        |   ...             |  ...     |  ...     |
+-------------+----------------------------+-------------------+----------+----------+

.. raw:: html

    <div style="text-align: justify;">
    Description of attributes for each column in the table file:
    <div style="margin-left: 40px;">
    <p></p>
    a. <span class="keyword-highlight">PDB</span>: provide a <span class="keyword-highlight">PDB</span> identifier that allow program can automatically download the <span class="keyword-highlight">PDB</span> file.
    <p></p>
    b. <span class="keyword-highlight">Amino Acid Substitution</span>: It consists of one-letter code of the wild-type amino acid, the sequential number of the mutation site, and the code of the mutant amino acid, for describing substitution of amino acids caused by the mutation.
    <p></p>
    d. <span class="keyword-highlight">Chain ID</span>: Indicate the protein chain where the mutation site is located.
    <p></p>
    e. <span class="keyword-highlight">pH</span>: Specify at which pH the mutation occurs. If you have no specific requirements or preferences regarding pH, you can simply specify it as 7.
    <p></p>
    f. <span class="keyword-highlight">T</span>: Specify at which temperature the mutation occurs. If you have no specific requirements or preferences regarding temperature, you can simply specify it as 25.
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
    This parameter requires user to provide a value of 0 or 1. The value of 0 means only predicting the ΔΔG for the mutations provided in the file, while the value of 1 means also predicting the ΔΔG for the reverse mutations of the mutations provided.
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
    This parameter requires user to provide an integer greater than 0 and less than 200. It represents the number of processes (multiprocessing) DDGWizard will use for calculating the optimal features.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>3. Output</h4>
    There will be an output csv file <span class="keyword-highlight">Pred_ddG.csv</span> located in <span class="keyword-highlight">DDGWizard/src/Pred_Res/</span>, which will record all prediction results.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>4. Notes</h4>
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
