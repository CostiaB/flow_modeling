def write_params(path, filename='start_params.py', **kwargs):
    line = '{} = {}\n'
    with open(path + '/' + filename, 'w') as f:
        for kw in kwargs:
            f.write(line.format(kw, kwargs[kw]))