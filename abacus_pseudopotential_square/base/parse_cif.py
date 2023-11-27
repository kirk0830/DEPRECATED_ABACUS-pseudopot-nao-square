"""Crsytal information file parser.

This module is used to parse the crystal information file (CIF) and return the
crystal information in dictionaries.

Methodology:
CIF files format varies from different software.
This module ONLY reads information starts everytime loop_ appears.

loop_ logic:
in loop_, if there is a line starts with "_", it means the title of the
information, if the next line also starts with "_", it means the information
will behind
"""
def parse(filename: str):

    with open(filename, 'r') as f:
        cif_contents = f.readlines()
    # concatenate the lines into a string
    concat_string = ""
    for line in cif_contents:
        if line.startswith("#"):
            continue # skip comments
        concat_string += line
    concat_string = concat_string.replace("\n", " ")
    # seperate the string into different loop_ sections
    loop_sections = concat_string.split("loop_")
    # parse each loop_ section
    cifs = []
    for section in loop_sections:
        cifs.append(parse_loop(section))
    return cifs

def parse_loop(concat_string: str) -> dict:
    """Parse the information after loop_.

    Args:
        concat_string: The string that contains the information in loop_.

    Returns:
        A dictionary that contains the information in loop_.
    """
    scattered_information = concat_string.split()
    if not scattered_information[0].startswith("_"):
        return # it is not a standard loop_ section, may be the super title
    # other normal case
    # 1. recombine contents in quotes into one item

    within_quotes = False
    content_in_quotes = ""
    recombined_information = []
    for word in scattered_information:
        if word.startswith("\'") or word.startswith("\"") or word.startswith(";"):
            within_quotes = True
            content_in_quotes += word.replace("\'", "").replace("\"", "").replace(";", "") + " "
            continue
        
        if word.endswith("\'") or word.endswith("\"") or word.endswith(";"):
            within_quotes = False
            content_in_quotes += word.replace("\'", "").replace("\"", "").replace(";", "")
            recombined_information.append(content_in_quotes)
            content_in_quotes = ""
            continue
            
        if within_quotes:
            content_in_quotes += " "+word
        else:
            recombined_information.append(word)
    print(recombined_information)
    # 2. identify cases between key-value and key-key-value-value
    attributes = [] # attributes of lines in order
    for word in recombined_information:
        if word.startswith("_"):
            attributes.append("key")
        else:
            attributes.append("value")
    print(attributes)
    # 3. seperate the key-value and key-key-value-value. assume key-value as a special case of key-key-value-value
    #    subloop is defined as a group of key-...-value
    subloops = {}
    index = 0
    cache_attribute = ""
    for ia, attribute in enumerate(attributes):
        if cache_attribute == attribute:
            if attribute == "key":
                subloops[index]["keys"].append(recombined_information[ia])
            else:
                subloops[index]["values"].append(recombined_information[ia])
        else:
            if attribute == "key": # case value-value-...-value-key
                index += 1
                subloops[index] = {
                    "keys": [recombined_information[ia]],
                    "values": [],
                }
            else: # case key-key-...-key-value
                subloops[index]["values"].append(recombined_information[ia])

            cache_attribute = attribute
    print(subloops)

    # 4. re-organize subloops
    cif = {}
    for index in subloops.keys():
        subloop = subloops[index] # get one subloop
        nkeys = len(subloop["keys"])
        key = ""
        print(subloop)
        for iv, value in enumerate(subloop["values"]):
            if int(iv/nkeys) == 0:
                key = subloop["keys"][iv%nkeys]
                cif[key] = []
            cif[key].append(value)
    print(cif)
    return cif

if __name__ == "__main__":

    print(parse("1534888.cif"))