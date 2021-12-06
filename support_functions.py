'''
Created on Jun 4, 2014

@author: levy
'''

import numpy as np


def path_summary(x):
    n = len(x)
    x1 = x[:(n-1)]
    x2 = x[1:]
    y = np.where(x1 != x2)[0]    
    start = [0] + list(y + 1)
    end = list(y) + [n-1]
    value = x[start]    
    return np.transpose(np.array((start, end, value)))

## loads a configuration file
## contains entries of the form:
## key=value
## and comment lines that are ignored.
## Expects string values in quotes, then will try to cast for int, then float.
def load_configuration_file(conf_filename, COMMENT_CHAR="#"):
    infile = file(conf_filename, 'r')
    
    ans = {}
    
    for line in infile:
        if line.startswith(COMMENT_CHAR):
            continue
        entry = [x.strip() for x in line.split("=")]
        if len(entry) != 2:
            continue
        ## if its a string with quotes, knock off the quotes
        if entry[1].startswith("'") or entry[1].startswith('"'):
            entry[1] = entry[1][1:-1]
        else:
            try:
                entry[1] = int(entry[1])
            except:
                try:
                    entry[1] = float(entry[1])
                except:
                    continue
        ans[entry[0]] = entry[1]
    
    return ans

def flatten_chararray(char_array):
    return char_array.view("S%d" % len(char_array))[0]

## nice function for printing a line
def tabprint(A, delim="\t"):
    return delim.join(map(str, A))

def tabprint_withformat(A, delim="\t", floatformat="%0.2f"):
    B = [floatformat % x if isinstance(x, float) else x for x in A]
    return delim.join(map(str, B))


def returnsValue(func, x):
    try:
        func(x)
        return True
    except:
        return False

def describe_array(ndarray, examples=1):
    for x, y in ndarray.dtype.descr:
        print tabprint([x, y, tabprint(ndarray[x][:examples], ",")])
        


def file2RA(filename, delimiter="\t", formats = None, headings=None):
    infile = file(filename, 'r')
    
    if headings==None:
        headings = infile.next().strip().split(delimiter)
    
    numCol = len(headings)
    
    isInt           = np.ones(numCol, dtype=bool)
    isFloat         = np.ones(numCol, dtype=bool)
    stringLength    = np.zeros(numCol, dtype=int)
    
    data = []
    for line in infile:
        entry = line.strip().split(delimiter)
        if len(entry) != numCol:
            print "skipping entry:",
            print line.strip()
        else:
            stringLength = [max(len(a), b) for a, b in zip(entry, stringLength)]
            isInt        = [(returnsValue(int, a) and b) for a,b in zip(entry, isInt)]
            isFloat        = [(returnsValue(float, a) and b) for a,b in zip(entry, isFloat)]
            data.append(entry)
    
    infile.close()
    
    if formats == None:
        formats = []
        for name, i, f, sl in zip(headings, isInt, isFloat, stringLength):
            if i:
                form = 'int'
            else:
                if f:
                    form = 'float'
                else:
                    form = 'S%i' % sl
            formats.append(form)
            
        for entry in data:
            for i in xrange(numCol):
                if formats[i] == 'int':
                    entry[i] = int(entry[i])
                elif formats[i] == 'float':
                    entry[i] = float(entry[i])
                    
    data2 = []
    for entry in data:
        data2.append(tuple(entry))
    
    DT = {'names':headings,
          'formats':formats}
    return np.array(data2, dtype=DT)


def files2RA(filename_list, delimiter="\t", formats = None):
    for filename in filename_list:
        infile = file(filename, 'r')    
        headings = infile.next().strip().split(delimiter)
        numCol = len(headings)
        
        isInt           = np.ones(numCol, dtype=bool)
        isFloat         = np.ones(numCol, dtype=bool)
        stringLength    = np.zeros(numCol, dtype=int)
        
        data = []
        for line in infile:
            entry = line.strip().split(delimiter)
            if len(entry) < numCol:
                print "skipping entry:",
                print line.strip()
            else:
                stringLength = [max(len(a), b) for a, b in zip(entry, stringLength)]
                isInt        = [(returnsValue(int, a) and b) for a,b in zip(entry, isInt)]
                isFloat        = [(returnsValue(float, a) and b) for a,b in zip(entry, isFloat)]
                data.append(entry)
        
        infile.close()
    
    if formats == None:
        formats = []
        for name, i, f, sl in zip(headings, isInt, isFloat, stringLength):
            if i:
                form = 'int'
            else:
                if f:
                    form = 'float'
                else:
                    form = 'S%i' % sl
            formats.append(form)
            
        for entry in data:
            for i in xrange(numCol):
                if formats[i] == 'int':
                    entry[i] = int(entry[i])
                elif formats[i] == 'float':
                    entry[i] = float(entry[i])
                    
    data2 = []
    for entry in data:
        data2.append(tuple(entry))
    
    DT = {'names':headings,
          'formats':formats}
    return np.array(data2, dtype=DT)
