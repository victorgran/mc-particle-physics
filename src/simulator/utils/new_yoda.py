# TODO: Open a yoda file without installing yoda.
# TODO: Also implement writing a YODA file here.


def writeYODA(filename: str):
    return


def readYODA(filename: str):
    with open(filename) as file:
        counter = 0
        for line in file:
            print(line)
            counter += 1
            if counter > 10:
                break
    return
