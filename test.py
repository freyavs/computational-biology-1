import re

from project3 import get_all_combinations as getallcombs3
from project1 import get_all_combinations as getallcombs1


def test_smallstring_sets_vs_tries():

    combs = getallcombs3("TESTBESTAND")
    print(len(combs))
    c = getallcombs1("TESTBESTAND")
    print(len(c))

    return combs.difference(c)



print(test_smallstring_sets_vs_tries)