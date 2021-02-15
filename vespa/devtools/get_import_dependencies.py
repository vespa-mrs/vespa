import importlib
import os
import pathlib
import re
import sys   

sys.setrecursionlimit(100000000)

dependenciesPaths = list()
dependenciesNames = list()
paths = sys.path


def main(path):
    """
    This recursively calls itself for whatever import modules it finds
    starting with 'path' and working its way down. It prints out the filename
    and module import paths that it finds/opens. 
    
    It returns a sorted list of module names that are the dependencies
    
    This is a pretty rough and ready implementation, but does a decent initial
    job of parsing the modules needed. 
    
    """
    try:
        print("Finding imports in '" + path + "':----------------------------------------------------------------------")

        file = open(path)
        contents = file.read()
        wordArray = re.split(" |\n", contents)

        currentList = list()
        nextPaths = list()
        skipWord = -1

        for wordNumb in range(len(wordArray)):
            word = wordArray[wordNumb]

            if wordNumb == skipWord:
                continue

            elif word == "from":
                item = wordArray[wordNumb + 1]
                if 'vespa.' in item:
                    currentList.append(item)
                skipWord = wordNumb + 2

            elif word == "import":
                item = wordArray[wordNumb + 1]
                if 'vespa.' in item:
                    currentList.append(item)

        currentList = set(currentList)
        for i in currentList:
            print(i)

        # print("Found imports in '" + path + "'")
        # print("Finding paths for imports in '" + path + "':")

        currentList2 = currentList.copy()
        currentList = list()

        for i in currentList2:
            if i in dependenciesNames:
                # print(i, "already found")
                pass

            else:
                dependenciesNames.append(i)

                try:
                    fileInfo = importlib.machinery.PathFinder().find_spec(i)
                    if fileInfo is None:
                        fileInfo = importlib.util.find_spec(i)
                        if fileInfo is None:
                            origin = 'None'
                        else:
                            origin = fileInfo.origin
                    else:
                        origin = fileInfo.origin

                    print(origin)
                    dependenciesPaths.append(origin)
                    currentList.append(origin)

                except AttributeError as e:
                    print("Hit Exception: AttributeError")
                    print(e)
                    print(i)
                    print(importlib.machinery.PathFinder().find_spec(i))
                    # print(red, "Odd noneType import called ", i, " in path ", path, end, sep='')


#        print("Found paths for imports in '" + path + "'")


        for fileInfo in currentList:
            main(fileInfo)

    except Exception as e:
        print(e)


if __name__ == "__main__":
    # args
    args = sys.argv
    print(args)

    if len(args) == 2:
        p = args[1]

    elif len(args) == 3:
        p = args[1]

        open(args[2], "a").close()
        sys.stdout = open(args[2], "w")

    else:
        print('Usage')
        print('PyDependencies <InputFile>')
        print('PyDependencies <InputFile> <OutputFile')

        sys.exit(2)

    if not os.path.exists(p):
        print("Path '" + p + "' is not a real path", sep='')
        #print(red, "Path '" + p + "' is not a real path", end, sep='')

    elif os.path.isdir(p):
        print("Path '" + p + "' is a directory, not a file", sep='')
        #print(red, "Path '" + p + "' is a directory, not a file", end, sep='')

    elif "".join(pathlib.Path(p).suffixes) != ".py":
        print("Path '" + p + "' is not a python file", sep='')
        #print(red, "Path '" + p + "' is not a python file", end, sep='')

    else:
        print("Path '" + p + "' is a valid python file", sep='')
        #print(green, "Path '" + p + "' is a valid python file", end, sep='')

        main(p)

    deps = set(dependenciesNames)

    deps = sorted(deps)

    print(deps)

    sys.exit()