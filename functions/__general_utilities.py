'''
this script includes some general utilities (not specific for a tool like McSnow, PAMTRA,...)
'''

import numpy as np
def gen_shortest_strings_without_duplicate(count_len=4):
    
    import itertools #used to generate all permutations

    #allocate list of strings
    if count_len==3: #for a string of 3 we have 238328 possible entries
        n_entries=238328
        count_str = ['___']*n_entries
    elif count_len==4: #for a string of 4 we have 14776336 possible entries
        n_entries=14776336
        count_str = ['___']*n_entries
    elif 14776336<number_ofSP:
        print "check number of SP (exceeding anticipated number in this adaption)"
        sys.exit(0)

    #generate list of strings with itertools
    #if __name__ == "__main__":
    chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

    i=0
    for item in itertools.product(chars, repeat=count_len):
        count_str[i] ="".join(item)
        i=i+1
    return count_str