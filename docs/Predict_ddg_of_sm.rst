.. _Predict ΔΔG for Saturation Mutagenesis:

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

Predict ΔΔG for Saturation Mutagenesis
=======================================

.. raw:: html

    <div style="text-align: justify;">
    This guide is intended to show users how to quickly use DDGWizard to predict ΔΔG for saturated mutagenesis.
    <p></p>
    Saturation mutagenesis represents mutating the original amino acid residue at the same mutation site to all possible amino acids. In practical applications, users often require predicting the ΔΔG of saturation mutagenesis at one or all amino acid sites, thereby assessing which mutations may enhance thermostability of the protein based on a wide range of possibilities.
    <p></p>
    To meet this practical user's need, we have prepared a program to help users quickly generate the needed <span class="keyword-highlight">csv</span> file for saturation mutagenesis. This file serves as input for DDGWizard to predict the ΔΔG of saturation mutagenesis.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>1. Running example</h4>
    <p>Similarly, we first provide two examples of running this program, followed by a detailed explanation of the program parameters. We selected the protein <span class="keyword-highlight">1SHG</span> as a case study.</p>
    <p></p>
    For predicting the ΔΔG of saturation mutagenesis at a single site (e.g. number 57 amino acid site), run the program with:
    <p></p>
    <div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ conda activate DDGWizard
    $ cd <b>&lt;/path/to/DDGWizard/&gt;</b>
    $ python Utility_Tool.py \
       --pdb_id 1SHG \
       --chain_id A \
       --site_number 57 \
       --wt_aa Y \
       --pH 7 \
       --T 25
    </pre>
    </div>
    </div>


.. raw:: html

    <div style="text-align: justify;">
    <p></p>
    For predicting the ΔΔG of full-site saturation mutagenesis, run the program with:
    <p></p>
    <div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ conda activate DDGWizard
    $ cd <b>&lt;/path/to/DDGWizard/&gt;</b>
    $ python Utility_Tool.py \
       --pdb_id 1SHG \
       --site_number all \
       --pH 7 \
       --T 25
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>2. Parameter details</h4>
    Below are the details of the parameters for the program of saturation mutagenesis:

.. raw:: html

    <div style="text-align: justify;">
    (1). <span class="keyword-highlight">--pdb_id</span>
    <p></p>
    Provide a <span class="keyword-highlight">PDB</span> identifier that allow program can automatically download the <span class="keyword-highlight">PDB</span> file.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    (2). <span class="keyword-highlight">--chain_id</span>
    <p></p>
    Indicate the protein chain where the mutation site is located.
    <p></p>
    If you intend to predict the ΔΔG of full-site saturation mutagenesis and the parameter <span class="keyword-highlight">--site_number</span> was provided with the value <span class="keyword-highlight">all</span>, you don't need to provide this parameter. The program will automatically match the chain identifier for all possible mutations.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    (3). <span class="keyword-highlight">--site_number</span>
    <p></p>
    This parameter indicates that you need to provide the site number of the predicted mutation.
    <p></p>
    If you intend to predict the ΔΔG of full-site saturation mutagenesis, please provide the value <span class="keyword-highlight">all</span>.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    (4). <span class="keyword-highlight">--wt_aa</span>
    <p></p>
    This parameter indicates that you need to provide the wild-type amino acid of the predicted mutation.
    <p></p>
    If you intend to predict the ΔΔG of full-site saturation mutagenesis and the parameter <span class="keyword-highlight">--site_number</span> was provided with the value <span class="keyword-highlight">all</span>, you don't need to provide this parameter. The program will automatically match the wild-type amino acid for all possible mutations.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    (5). <span class="keyword-highlight">--pH</span>
    <p></p>
    This parameter indicates that you need to specify at which pH you want to predict the ΔΔG for the mutations.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    (6). <span class="keyword-highlight">--T</span>
    <p></p>
    This parameter indicates that you need to specify at which temperature you want to predict the ΔΔG for the mutations.
    <p></p>
    <div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>3. Output</h4>
    The program will generate an output csv file <span class="keyword-highlight">Pred.csv</span> located in <span class="keyword-highlight">DDGWizard/src/</span>.
    <p></p>
    This <span class="keyword-highlight">csv</span> file can be directly used as input for the DDGWizard prediction program, enabling quick preparation for ΔΔG prediction of saturation mutagenesis:
    <p></p>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ conda activate DDGWizard
    $ cd <b>&lt;/path/to/DDGWizard/&gt;</b>
    $ python Predict_ddG_Executable.py \
        --pred_dataset_path ./src/Pred.csv \
        --db_folder_path <b>&lt;/folder/to/save/Blast_database/&gt;</b> \
        --db_name <b>&lt;the_name_to_assign_for_Blast database&gt;</b> \
        --if_reversed_data 0 \
        --blast_process_num 4 \
        --mode whole \
        --process_num 4
    </pre>
    </div>
    </div>

