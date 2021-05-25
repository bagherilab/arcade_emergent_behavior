from .utilities import load, load_tar
import tarfile

class PARAMETER_SENSITIVITY():
    NAME = "PARAMETER_SENSITIVITY"

    PARAMETERS = ['max_height','meta_pref','migra_threshold']

    PERCENTS = ['000','010','020','030','040','050','060','070','080','090','100',
        '110','120','130','140','150','160','170','180','190','200']

    @staticmethod
    def run(output_path, input_path, func, name=NAME,
            parameters=PARAMETERS, percents=PERCENTS, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for param in parameters:
            for perc in percents:
                code = f"_{param}_{perc}"
                infile = f"{input_path}{name}/{name}_{param}_{perc}.pkl"
                print(f"{name} : {code}")

                loaded = load(infile)
                func(*loaded, outfile, code, exclude=[-1], timepoints=timepoints, seeds=seeds)

    @staticmethod
    def loop(output_path, func1, func2, extension, name=NAME,
             parameters=PARAMETERS, percents=PERCENTS, timepoints=[]):
        outfile = f"{output_path}{name}/{name}"
        out = { "data": [] }
        tar = load_tar(outfile, extension)

        for t in timepoints:
            for param in parameters:
                for perc in percents:
                    code = f"_{param}_{perc}"
                    func1(outfile, out, { "time": t, "param": param, "perc": perc }, extension, code, tar=tar)

        func2(outfile, extension, out)

    @staticmethod
    def load(output_path, input_path, func, extension="", name=NAME,
             parameters=PARAMETERS, percents=PERCENTS, timepoints=[], seeds=[]):
        outfile = f"{output_path}{name}/{name}"

        for param in parameters:
            for perc in percents:
                code = f"_{param}_{perc}"
                infile = f"{input_path}{name}{extension}/{name}_{param}_{perc}{extension}.tar.xz"
                print(f"{name} : {code}")

                tar = tarfile.open(infile)
                func(tar, timepoints, { "param": param, "perc": perc }, outfile, code)
