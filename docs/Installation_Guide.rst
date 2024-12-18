.. _installation_guide:

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

Installation Guide
==================

.. raw:: html

    <div style="text-align: justify;">
    DDGWizard consists of 3 components: the feature calculation pipeline, that processes raw ΔΔG data and outputs feature-enriched ΔΔG data with 1547 features; the DDGWizard dataset, including 15752 ΔΔG data; and the accurate ΔΔG prediction model.
    <p></p>
    This section explains how to install dependencies for using the DDGWizard's application (there is no need to install anything to access the DDGWizard dataset; it can be directly downloaded).
    <p></p>
    <h4>Installation prerequisites:</h4>
    CentOS 7 or Ubuntu system; GCC version higher than 4.8.5; Conda version higher than 23.0; Git.
    <p></p>
    </div>

.. _`the Characterization part`:

Feature Calculation Pipeline (for Generating Feature-Enriched ΔΔG Data)
------------------------------------------------------------------------------

.. raw:: html

    <div style="text-align: justify;">
    This subsection is for users who need to use the feature calculation pipeline. It can assist users in processing input raw ΔΔG data and outputting feature-enriched new data, including 1574 features that completed calculations. It can facilitate further analysis, feature selection, and machine learning.
    <p></p>
    The installation steps are as follows.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>1. Git clone the DDGWizard repository</h4>
    <p></p>
    </div>

.. code-block::

    $ git clone https://github.com/bioinfbrad/DDGWizard.git

.. raw:: html

    <div style="text-align: justify;">
    <h4>2. Config and install conda virtual environment</h4>
    <p></p>
    <p>There is an <span class="keyword-highlight">Environment.yml</span> file located in the path <span class="keyword-highlight">DDGWizard/src</span>, which is the Conda environment configuration file.</p>
    <p></p>
    <p>Open this file with your text editor (e.g., nano, vim, vi, etc.). Here we use vi as an example:</p>
    <p></p>
    </div>

.. code-block::

    $ cd DDGWizard/src/
    $ vi Environment.yml

.. raw:: html

    <div style="text-align: justify;">
    <p>Modify the <span class="keyword-highlight">prefix</span>, <b>which is on the last line</b>. <b>Change the prefix to your local <span class="keyword-highlight">conda envs folder</span>.</b></p>
    <p></p>
    After changing, the <span class="keyword-highlight">prefix</span> should be <span class="keyword-highlight">prefix: <b>&lt;the path to your conda envs folder&gt;</b>/DDGWizard</span>.
    <p></p>
    If you don't know how to find the path to local <span class="keyword-highlight">conda envs folder</span>, you can use command:
    <p></p>
    </div>

.. code-block::

     $ conda info --envs

.. raw:: html

    <div style="text-align: justify;">
    <p>Once user have changed the <span class="keyword-highlight">prefix</span> of <span class="keyword-highlight">Environment.yml</span> file, please use Conda commands to create a Conda virtual environment and install dependencies. This may take some time.</p>
    <p></p>
    </div>

.. code-block::

     $ conda env create -f Environment.yml

.. raw:: html

    <div style="text-align: justify;">
    <h4>3. Download NCBI-BLAST-2.13.0+</h4>
    <p></p>
    Users need to download the NCBI-BLAST-2.13.0 program for allowing DDGWizard to carry out multiple sequence alignment (MSA). Please visit <a href="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/">Download NCBI-BLAST-2.13.0+</a> to download the <span class="keyword-highlight">ncbi-blast-2.13.0+-x64-linux.tar.gz</span> file. Copy this compressed file to the path <span class="keyword-highlight">DDGWizard/bin/ncbi_blast_2_13_0+/</span> and extract it. Use the following commands:
    <p></p>
    </div>

.. code-block::

     $ wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
     $ cp ncbi-blast-2.13.0+-x64-linux.tar.gz DDGWizard/bin/ncbi_blast_2_13_0+/
     $ cd DDGWizard/bin/ncbi_blast_2_13_0+/
     $ tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz
     $ cp -r ncbi-blast-2.13.0+/* .

.. raw:: html

    <div style="text-align: justify;">
    NCBI-BLAST-2.13.0+ is a "United States Government Work" under the terms of the United States Copyright Act. Please read and accept the license file in its folder before proceeding further.
    <p></p>
    </div>

.. raw:: html

   <div style="text-align: justify;">
   <h4>4. Configure Modeller</h4>
   <p></p>
   The Modeller software is used for homology or comparative modeling of protein three-dimensional structures. In DDGWizard, Modeller is used to construct PDB protein structure files of mutations based on the user's input of wild-type PDB protein structure files.
   <p></p>
   Modeller has already been installed when creating Conda environment. But to allow our program to call it, you need to have a license of the Modeller and configure it.
   <p></p>
   Please enter <a href="https://salilab.org/modeller/registration.html">Official Modeller Website</a> to register an account. Modeller use "Academic End-User Software License Agreement for MODELLER" terms. Please follow their instructions, read and accept the terms to obtain a license.
   <p></p>
   Then input the license into installed Modeller's configuration file. You can find it under the <span class="keyword-highlight">Conda envs folder</span>.
   <p></p>
   Enter your local <span class="keyword-highlight">Conda envs folder</span>, and open the Modeller's configuration file:
   <p></p>
   </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cd <b>&lt;the path to your conda envs folder&gt;</b>
    $ vi DDGWizard/lib/modeller-10.6/modlib/modeller/config.py
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <p>Replace the XXXX to your license. Save and close it.</p>
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4> To use DDGWizard feature calculation pipeline, the following software is optional (step 5-11) and not required to be installed (if certain software is not installed, the feature values it calculates will not be output).</h4>
    <p></p>
    <h4>Before running, please don't forget to make sure the programs of the DDGWizard have the executable permission (step 12). Return to the DDGWizard program folder and execute the command:</h4>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cd DDGWizard/
    $ chmod -R +x .
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4> (Optional) 5. Download FoldX 5.0</h4>
    <p></p>
    Users can download the FoldX 5.0 program for allowing DDGWizard to calculate energy terms of proteins. FoldX has academic version and commercial version. To use it in DDGWizard, academic version is enough. Please visit <a href="https://foldxsuite.crg.eu/academic-license-info">Apply for FoldX 5.0</a> to register an account, read and accept "FoldX Academic License" terms to download the <span class="keyword-highlight">foldx5Linux64.zip</span> file. Copy this compressed file to the path <span class="keyword-highlight">DDGWizard/bin/FoldX_5.0/</span> and extract it. Use the following commands:
    <p></p>
    </div>

.. code-block::

     $ cp foldx5Linux64.zip DDGWizard/bin/FoldX_5.0/
     $ cd DDGWizard/bin/FoldX_5.0/
     $ unzip foldx5Linux64.zip

.. raw:: html

    <div style="text-align: justify;">
    <h4> (Optional) 6. Download Ring 3.0</h4>
    <p></p>
    Users can download the Ring 3.0 application for allowing DDGWizard to calculate residue interaction information. Please visit <a href="https://biocomputingup.it/services/download/">Apply for Ring 3.0</a> to apply. Please read and accept the license of Ring 3.0 to obtain the <span class="keyword-highlight">ring-3.0.0.tgz</span> file. Copy this compressed file to the path <span class="keyword-highlight">DDGWizard/bin/ring-3.0.0/</span> and extract it. Use the following commands:
    <p></p>
    </div>

.. code-block::

     $ cp ring-3.0.0.tgz DDGWizard/bin/ring-3.0.0/
     $ cd DDGWizard/bin/ring-3.0.0/
     $ tar -zxvf ring-3.0.0.tgz
     $ cp -r ./ring-3.0.0/* .

.. raw:: html

    <div style="text-align: justify;">
    <h4> (Optional) 7. Download DisEMBL</h4>
    <p></p>
    Users can download the DisEMBL program for allowing DDGWizard to count disorder information of proteins. Please visit <a href="https://zenodo.org/records/14246673">Download the DisEMBL</a> to download the <span class="keyword-highlight">DisEMBL-1.4.tgz</span> file. Copy this compressed file to the path <span class="keyword-highlight">DDGWizard/bin/DisEMBL_1_4/</span> and extract it. Use the following commands:
    <p></p>
    </div>

.. code-block::

     $ cp DisEMBL-1.4.tgz DDGWizard/bin/DisEMBL_1_4/
     $ cd DDGWizard/bin/DisEMBL_1_4/
     $ tar -zxvf DisEMBL-1.4.tgz
     $ cp -r ./DisEMBL-1.4/* .

.. raw:: html

    <div style="text-align: justify;">
    DisEMBL uses GPL 2.0 open-source license. Please read and accept the license file in its folder before proceeding further.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4> (Optional) 8. Configure DSSP</h4>
    The DSSP is used to calculate the RSA (relative surface area) and secondary stuctures of <span class="keyword-highlight">PDB</span> files.
    <p></p>
    To allow DDGWizard use DSSP, please enter your local <span class="keyword-highlight">Conda envs folder</span>, then enter <span class="keyword-highlight">bin folder</span>, and copy <span class="keyword-highlight">mkdssp</span> as <span class="keyword-highlight">dssp</span>:
    <p></p>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cd <b>&lt;the path to your conda envs folder&gt;</b>
    $ cd DDGWizard/bin/
    $ cp mkdssp dssp
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>(Optional) 9. Install Bio3D</h4>
    Users can install the Bio3D package for allowing DDGWizard to calculate atomic fluctuations based on NMA (normal mode analysis). It requires users have <span class="keyword-highlight">R</span> as prerequisites (it can be downloaded and installed from <a href="https://cran.r-project.org/">Official R Website</a>). Then please use following commands to install Bio3d package:
    <p></p>
    </div>

.. code-block::

    $ R
    install.packages("bio3d")

.. raw:: html

    <div style="text-align: justify;">
    <h4>(Optional) 10. Download PROFbval</h4>
    PROFbval relies on the Ubuntu environment. To address cross-platform compatibility, we have created container images for easy download by users. This requires users have Docker or Singularity as a prerequisite.
    <p></p>
    Please download the following two files: <span class="keyword-highlight">myprof.tar</span> (128MB) and <span class="keyword-highlight">myprof.sif</span> (360MB) from <a href="https://zenodo.org/records/12817843">https://zenodo.org/records/12817843</a>, and copy them to the path: <span class="keyword-highlight">DDGWizard/src/Prof_Source</span>:
    <p></p>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cp <b>&lt;the path to myprof.tar&gt;</b>/myprof.tar DDGWizard/src/Prof_Source
    $ cp <b>&lt;the path to myprof.sif&gt;</b>/myprof.sif DDGWizard/src/Prof_Source
    </pre>
    </div>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    DDGWizard will automatically call the programs within the container images. Users only need to have either Docker or Singularity. If the user chooses Docker, an additional step is required:
    <p></p>
    </div>

.. code-block::

    $ docker load -i DDGWizard/src/Prof_Source/myprof.tar

.. raw:: html

    <div style="text-align: justify;">
    PROFbval uses GPL 3.0+ open-source license. Please read and accept its license before proceeding further.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4> (Optional) 11. Download SIFT 6.2.1</h4>
    <p></p>
    Users can download the SIFT 6.2.1 program for allowing DDGWizard to predict impact of amino acid substitution on protein function. Please visit <a href="https://sift.bii.a-star.edu.sg/www/code.html">Download SIFT 6.2.1</a> to download the <span class="keyword-highlight">sift6.2.1.tar.gz</span> file. Copy this compressed file to the path <span class="keyword-highlight">DDGWizard/bin/sift6_2_1/</span> and extract it. Use the following commands:
    <p></p>
    </div>

.. code-block::

     $ wget https://s3.amazonaws.com/sift-public/nsSNV/sift6.2.1.tar.gz
     $ cp sift6.2.1.tar.gz DDGWizard/bin/sift6_2_1/
     $ cd DDGWizard/bin/sift6_2_1/
     $ tar -zxvf sift6.2.1.tar.gz
     $ cp -r sift6.2.1/* .

.. raw:: html

    <div style="text-align: justify;">
    SIFT 6.2.1 uses non-commercial license. Please read and accept the license file in its folder before proceeding further.
    <p></p>
    </div>

.. raw:: html

    <div style="text-align: justify;">
    <h4>12. Make sure the programs of the DDGWizard have the executable permission</h4>
    The programs of DDGWizard need the executable permission to run.
    <p></p>
    Return to the DDGWizard program folder and execute the command:
    <p></p>
    </div>

.. raw:: html

    <div class="highlight-default notranslate">
    <div class="highlight">
    <pre style="overflow: scroll">
    $ cd DDGWizard/
    $ chmod -R +x .
    </pre>
    </div>
    </div>

.. _`the Prediction Part`:

ΔΔG Prediction Model (for Predicting ΔΔG)
-----------------------------------------------

.. raw:: html

    <div style="text-align: justify;">
    This subsection is for users who need to use the ΔΔG prediction model.
    <p></p>
    To use DDGWizard's ΔΔG prediction model, users are required to complete steps 1-8 (these are no longer optional) and execute step 12 of Feature Calculation Pipeline's installation part. Steps 9-11 are not required.
    <p></p>
    </div>

