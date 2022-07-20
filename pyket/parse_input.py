#inspired by https://github.com/patrickmineault/zipf/blob/master/zipf/parse_text.py
import pathlib
import pandas as pd 

def read_photon_bath(f, clean_text=False):
    """
    Parses the photon bath input file.

    Arguments:
        f: an open file handle
        clean_text (optional): a Boolean, if true, filters comments starting with #
    
    Returns:
        A dict, keyed by photon frequency, with photon coupling as value. All items are strings
    """
    text = f.read()
    if clean_text:
        text = _clean_photon_bath(text)
    chunks = [pair.partition(':') for pair in text.split() if pair]
    photon_dict = {pair[0]:pair[2] for pair in chunks if not pair == ('','',':')}
    return photon_dict


def _clean_photon_bath(text):
    """
    Recursivly find lines starting with '#' and removes them. Also removes empty lines.
    """
    comment_start = text.find('#')
    #base case
    if(comment_start < 0 ):
        return text.strip()
    else:
        comment_end = text.find('\n',comment_start)
        text = text[:comment_start] + text[comment_end:]
        return _clean_photon_bath(text)
    

    