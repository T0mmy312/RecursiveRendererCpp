
#
# this just builds the finished project, that everything is in one header file just so its simpler.
#

toInclude = []

def getContent(filename):
    ret = ""
    with open(filename, "r") as file:
        con = file.read().split("\n")
        for line in con:
            if '#include "' in line:
                ret += getContent("include/" + line[10:len(line) - 1:1]) + "\n"
            elif "#pragma once" in line:
                continue
            elif "#include <" in line:
                found = False
                for include in toInclude:
                    if line in include:
                        found = True
                if not found:
                    toInclude.append(line + "\n")
            else:
                ret += line + "\n"
        file.close()
    return ret

content = getContent("include/pngWrite.h") + "\n#endif"
content = "\n#pragma once\n\n#ifndef _COMPLETERENDERER_H_\n#define _COMPLETERENDERER_H_\n\n" + "".join(toInclude) + content
with open("completeRenderer.h", "w") as file:
    file.write(content)
    file.close()