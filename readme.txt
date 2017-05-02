# readme.txt

sdssdisp is created by Ai-Lei Sun in August 2015 to do broadband image operations on SDSS spec/photo objects. The main purpose is to study the emission line region properties of luminous Type 2 AGN. 

Example Usage: 

    To make SDSS/data/magellan/ and the contents therein: 

        > python main_magellan.py

    To make SDSS/data/mullaney/ and the contents therein: 

        > python main_mullaney.py

    To make SDSS/data/SDSS/SDSSJ1356+1026/ and the contents therein: 

        > python main_1356.py


File descriptions:


    main scripts:

        main.py             (under construction)
        main_1356.py        
        main_magellan.py    
        main_mullaney.py    

    class definition:

        class_obsobj.py     Class of SDSS Object identified by ra, dec 


    functions:

        alignstamp.py       make aligned stamp multi-band images of an 
                            SDSS object

        imagedisp_util.py   A collection of functions for image displays, 
                            e.g. objw_HumVIgriimages() that makes .pgn 
                            color images.

        subtractimg.py      To produce continuum subtracted broadband images 
                            and [OIII] intensity maps. It uses fromspec 
                            functions to determine the raito. 

        fromspec.py         Download and store SDSS spectrum of objects, 
                            calculate the continuum flux level at different 
                            bands, and return the ratio of continuum flux 
                            level between band1 and band2. 

        measurenebula.py    measure the nebula properties, like size, area, 
                            luminosities from SDSS images, and store them in 
                            the measure....txt file under each batch 
                            directory.  

        run_batch.py        do all the above image operations and 
                            measurements on a batch of SDSS objects. 


    special operation: 

        magellan_getslit.py put artificial slits on the SDSS images that 
                            represent Magellan observations, extract broad-
                            band-constructed line intensities along those 
                            slits, and compare with the Magellan 
                            observations.   


