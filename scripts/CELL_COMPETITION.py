from .utilities import load, load_tar
import tarfile

class CELL_COMPETITION():
    NAME = "CELL_COMPETITION"

    PARAMETERS = ['max_height','meta_pref','migra_threshold']

    PERCENTS = ['050','060','070','080','090','100','110','120','130','140','150']

    INITIALIZATIONS = [('000','100'), ('010','090'), ('020','080'),
        ('030','070'), ('040','060'), ('050','050'), ('060','040'),
        ('070','030'), ('080','020'), ('090','010'), ('100','000')]

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            parameters=PARAMETERS, percents=PERCENTS, initializations=INITIALIZATIONS,
            timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for param in parameters:
            for perc in percents:
                for init_basal, init_modified in initializations:
                    code = f"_{param}_{perc}_{init_basal}_{init_modified}"
                    infile = f"{input_path}{name}/{name}_{param}_{perc}_{init_basal}_{init_modified}.pkl"
                    print(f"{name} : {code}")

                    loaded = load(infile)
                    func(*loaded, outfile, code, exclude=[-1], timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             parameters=PARAMETERS, percents=PERCENTS, initializations=INITIALIZATIONS,
             timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for param in parameters:
                for perc in percents:
                    for init_basal, init_modified in initializations:
                        code = f"_{param}_{perc}_{init_basal}_{init_modified}"
                        func1(outfile, out, { "time": t, "param": param, "perc": perc, "init": f"{init_basal}-{init_modified}" }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             parameters=PARAMETERS, percents=PERCENTS, initializations=INITIALIZATIONS,
             timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for param in parameters:
            for perc in percents:
                for init_basal, init_modified in initializations:
                    code = f"_{param}_{perc}_{init_basal}_{init_modified}"
                    infile = f"{input_path}{name}{extension}/{name}_{param}_{perc}_{init_basal}_{init_modified}{extension}.tar.xz"
                    print(f"{name} : {code}")

                    tar = tarfile.open(infile)
                    func(tar, timepoints, { "param": param, "perc": perc, "init": f"{init_basal}-{init_modified}" }, outfile, code)
