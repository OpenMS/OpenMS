import sys
import os.path
import pyopenms as pms
import pprint

def addDataProcessing(obj, params, action):
    if isinstance(obj, pms.MSExperiment):
        result = pms.MSExperiment()
        for spec in obj:
            spec = _addDataProcessing(spec, params, action)
            result.push_back(spec)
    else:
        result = _addDataProcessing(obj, params, action)
    return result

def _addDataProcessing(item, params, action):
    dp = item.getDataProcessing()
    p = pms.DataProcessing()
    p.setProcessingActions(set([action]))
    sw = p.getSoftware()
    sw.setName(os.path.basename(sys.argv[0]))
    sw.setVersion(pms.VersionInfo.getVersion())
    p.setSoftware(sw)
    p.setCompletionTime(pms.DateTime.now())

    for k, v in params.asDict().items():
        p.setMetaValue("parameter: "+k, v)

    dp.append(p)
    item.setDataProcessing(dp)
    return item

def writeParamsIfRequested(args, params):
    if args.write_dict_ini or args.write_ini:
        if args.write_dict_ini:
            with open(args.write_dict_ini, "w") as fp:
                pprint.pprint(params.asDict(), stream=fp)
        if args.write_ini:
            fh = pms.ParamXMLFile()
            fh.store(args.write_ini, params)
        return True
    return False

def updateDefaults(args, defaults):
    if args.ini:
        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(args.ini, param)
        defaults.update(param)
    elif args.dict_ini:
        with open(args.dict_ini, "r") as fp:
            try:
                dd = eval(fp.read())
            except:
                raise Exception("could not parse %s" % args.dict_ini)
        defaults.updateFrom(dd)
