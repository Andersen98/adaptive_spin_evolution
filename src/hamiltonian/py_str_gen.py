from string import Template
if __name__=="__main__":
    out = ""
    s = Template('%|${mode}$$-d|:%|${level}$$-d|%|${space}t|')
    for i in range(10):
        d = dict(mode=str(2*i+1),level=str(2*i+2),space=str(8*(i+1)))
        out = out + s.substitute(d)
    out = out + "\\n"
    print(out)

    out = ""
    
    s = Template('fmtr % mode_level_array[${idx}].first;\nfmtr % mode_level_array[${idx}].second;\n')
    for i in range(10):
        d = dict(idx ="i*modes_per_line +"+ str(i))
        out = out + s.substitute(d)
    print(out)
