Example notebooks
=================
We have a set of notebooks that explain and show the usage of the SSPC features. The notebooks are also stored in the `examples directory in the repository <https://gitlab.com/dhendriks/syntheticstellarpopconvolve/-/tree/master/examples>`_

.. raw:: html


    <style>
        .flip-card {
            background-color: transparent;
            width: 220px;
            height: 300px;
            perspective: 1000px;
            margin: 10px;
        }

        .flip-card-inner {
            position: relative;
            width: 100%;
            height: 100%;
            text-align: center;
            transition: transform 0.8s;
            transform-style: preserve-3d;
        }

        .flip-card:hover .flip-card-inner {
            transform: rotateY(180deg);
        }

        .flip-card-front, .flip-card-back {
            position: absolute;
            width: 100%;
            height: 100%;
            backface-visibility: hidden;
            border: 2px solid #ccc;
            border-radius: 10px;
            overflow: hidden;
        }

        .flip-card-front img,
        .flip-card-back img {
            width: 100%;
            height: 100%;
            object-fit: contain;
        }

        .flip-card-back {
            transform: rotateY(180deg);
        }

        .flip-gallery {
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            justify-content: center;
        }

        .flip-caption {
            margin-top: 5px;
            font-size: 14px;
            text-align: center;
            color: #333;
        }
    </style>

    <div class="flip-gallery">

        <div class="flip-card">
            <a href="examples/notebook_example_GW_merger_rate_density.html">
                <div class="flip-card-inner">
                    <div class="flip-card-front">
                        <img src="_static/notebook_cover_images/GW_notebook_cover_image_front.png" alt="Front Cover">
                    </div>
                    <div class="flip-card-back">
                        <img src="_static/notebook_cover_images/GW_notebook_cover_image_back.png" alt="Back Cover">
                    </div>
                </div>
            </a>
            <div class="flip-caption">GW merger rate notebook</div>
        </div>

        <div class="flip-card">
            <a href="examples/notebook_example_LISA_UCB.html">
                <div class="flip-card-inner">
                    <div class="flip-card-front">
                        <img src="_static/notebook_cover_images/LISA_cover_image_front.png" alt="Front Cover">
                    </div>
                    <div class="flip-card-back">
                        <img src="_static/notebook_cover_images/LISA_cover_image_front.png" alt="Back Cover">
                    </div>
                </div>
            </a>
            <div class="flip-caption">LISA double white-dwarf population notebook</div>
        </div>

    </div>











.. toctree::
    :maxdepth: 2
    :caption: Contents:

    examples/notebook_convolution_tutorial.ipynb
    examples/notebook_convolution_star_formation_functions.ipynb
    examples/notebook_convolution_use_cases.ipynb
    examples/Background.ipynb

    examples/notebook_tutorial_persistent_data_and_previous_convolution_results.ipynb
    examples/notebook_tutorial_convolution_by_sampling.ipynb
    examples/notebook_tutorial_convolution_by_integration.ipynb
    examples/notebook_tutorial_convolution_on_the_fly.ipynb
    examples/notebook_tutorial_inflate_ensemble.ipynb
    examples/notebook_tutorial_convolve_binned_data.ipynb

    examples/notebook_example_bincodex.ipynb
    examples/notebook_example_GCE.ipynb
    examples/notebook_example_GW_merger_rate_density.ipynb
    examples/notebook_example_LISA_UCB.ipynb
    examples/notebook_example_orbit_integration_and_supernova_kicks.ipynb
    examples/notebook_example_GAIA_populations.ipynb
