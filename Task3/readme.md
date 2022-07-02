# TASK 3: Assigned 26062022

    Redo both sets of experiments (Task 1 and Task 2) after splitting the data energy wise (0.3-1 TeV, 1-3 TeV, 3-10 TeV)
    After that for method 2 (when counting no of matches within the neutrino error circle) instead of counting the background events by looking at the no of neutrino-MSP pairs within the same cos(error region) offset by 5 degrees,
    Follow the same method as in        https://arxiv.org/pdf/2112.11375v1.pdf      which is briefly summarized below
        1. Generate a synthetic catalog of MSPs uniformly distributed in the same RA and DEC range as observed population
        2, Count no of coincidences  within neutrino error circle for the synthetic catalog in (1). 
        3. Repeat step (1) and (2) 100 times and take the average.


## RESPONSE

    Used multithreading for this code as my code took too long in serial execution.
    Without multithreading, I ran the code yesterday for around 6 hrs without getting an output.
    With multithreading, I got the output in 131 mins(I used 12 Threads: One thread each for the 3 methods: 2b, 2c, 2d, so total 3 threads, and 
                                                                                each iteration calculates results for 4 synthetic sets of coordinates so $3*4=12$ threads!
    To save memory, I stored all the required methods used in method/task 2 in a separate .py file called load.py,  and called functions from it when required.
    To further decrease the execution time, I did the following:             
                    I did not plot the cosine(space angle) distribution as I did in method/task 2b.
                    I did not store and generate a new catalogue of matches as I did in method/tasks 2c and 2d but
                        just found the no.of matches satisfying their respective conditions.
