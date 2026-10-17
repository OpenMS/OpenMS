.. _workflows:

Workflow Editor
===============

Which workflow environment to choose for running OpenMS tools?

You can run OpenMS TOPP tools from the command line using your custom scripts, or use powerful
workflow systems designed to make workflow creation and maintenance more fun:


.. grid:: 2
    :gutter: 2

    .. grid-item-card:: Nextflow
        :img-top: /_images/introduction/nextflow2014_no-bg.png
        :link: /getting-started/nextflow-get-started
        :link-type: doc
        :class-card: only-light
        :class-img-top: halfwidth
        :columns: 12

        Script/DSL-based workflow language, executor and utilities (such as the browser based launcher and supervisor nf-tower). Automatically runs on various different cloud (AWS, Google, ...) and HPC environments (SLURM, LFS, Kubernetes, ...). It is recommended to use our ready-made nf-core compatible workflows for ease of use via the browser-based configuration and launcher.

    .. grid-item-card:: Nextflow
        :img-top: /_images/introduction/nextflow2014_no-bg-bright.png
        :link: /getting-started/nextflow-get-started
        :link-type: doc
        :class-card: only-dark
        :class-img-top: halfwidth
        :columns: 12

        Script/DSL-based workflow language, executor and utilities (such as the browser based launcher and supervisor nf-tower). Automatically runs on various different cloud (AWS, Google, ...) and HPC environments (SLURM, LFS, Kubernetes, ...). It is recommended to use our ready-made nf-core compatible workflows for ease of use via the browser-based configuration and launcher.


    .. grid-item-card:: Galaxy
        :img-top: /_images/introduction/galaxy_project_logo.png
        :link: /getting-started/galaxy-get-started
        :link-type: doc
        :class-card: only-light
        :class-img-top: halfwidth
        :columns: 12

        Server and browser-based interactive workflow editor and runner. A public server instance can be used for testing and smaller experiments. Provides nice guided tutorials.

    .. grid-item-card:: Galaxy
        :img-top: /_images/introduction/galaxy_project_logo_white.png
        :link: /getting-started/galaxy-get-started
        :link-type: doc
        :class-card: only-dark
        :class-img-top: halfwidth
        :columns: 12

        Server and browser-based interactive workflow editor and runner. A public server instance can be used for testing and smaller experiments. Provides nice guided tutorials.


    .. grid-item-card:: TOPPAS
        :img-top: /_images/introduction/TOPPAS_logo_white.png
        :link: /getting-started/toppas-get-started
        :link-type: doc
        :class-card: only-light
        :class-img-top: halfwidth
        :columns: 12

        OpenMS' build-in workflow system, with limited capabilities but easy to use and tailored to TOPP tools.

    .. grid-item-card:: TOPPAS
        :img-top: /_images/introduction/TOPPAS_logo_dark.png
        :link: /getting-started/toppas-get-started
        :link-type: doc
        :class-card: only-dark
        :class-img-top: halfwidth
        :columns: 12

        OpenMS' build-in workflow system, with limited capabilities but easy to use and tailored to TOPP tools.

    .. grid-item-card:: CWL
        :img-top: /_images/introduction/CWL-4K.png
        :link: /getting-started/cwl-get-started
        :link-type: doc
        :class-card: only-light
        :class-img-top: halfwidth
        :columns: 12

        OpenMS tools can self-generate CWL tool descriptions. As of OpenMS 3.6, TDL support is disabled by default; official binary releases include it. Rebuild with ``-DENABLE_TDL=ON`` to enable CWL generation.

    .. grid-item-card:: CWL
        :img-top: /_images/introduction/CWL-4K.png
        :link: /getting-started/cwl-get-started
        :link-type: doc
        :class-card: only-dark
        :class-img-top: halfwidth
        :columns: 12

        OpenMS tools can self-generate CWL tool descriptions. As of OpenMS 3.6, TDL support is disabled by default; official binary releases include it. Rebuild with ``-DENABLE_TDL=ON`` to enable CWL generation.

.. toctree::
    :maxdepth: 1
    :hidden:

    /getting-started/nextflow-get-started.md
    /getting-started/galaxy-get-started.md
    /getting-started/toppas-get-started.md
    /getting-started/cwl-get-started.md
