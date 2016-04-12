## Model that holds the gene and genome data
##  for the gene search

GeneSearchModel =
    vfs: getEndpoint(url="data/vf")
    amrs: getEndpoint(url="data/amr")
    genomes: @genomedata = getEndpoint(url="data/meta")

    vfList: {}
    amrList: {}
    genomeList: []
    selectedVF: []
    selectedAMR: []
    selectedGenomes: []

    setLists: ->
        vfdata = GeneSearchModel.vfs()
        amrdata = GeneSearchModel.amrs()
        GeneSearchModel.vfList["headers"] = vfdata.headers
        GeneSearchModel.amrList["headers"] = amrdata.headers
        for gene in vfdata.rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
            GeneSearchModel.vfList["rows"] = vfdata.rows
        for gene in amrdata.rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
            GeneSearchModel.amrList["rows"] = amrdata.rows

GeneSearchModel.setLists()