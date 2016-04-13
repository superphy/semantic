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

    getSelectedVF: ->
        GeneSearchModel.selectedVF = []
        for row in GeneSearchModel.vfList.rows when row.selected()
            GeneSearchModel.selectedVF.push(row.Gene_Name)


    getSelectedAMR: ->
        GeneSearchModel.selectedAMR = []
        for row in GeneSearchModel.amrList.rows when row.selected()
            GeneSearchModel.selectedAMR.push(row.Gene_Name)

    reset: ->
        console.log("Resetting...")
        GeneSearchModel.selectedVF = []
        GeneSearchModel.selectedAMR = []
        for gene in GeneSearchModel.vfList.rows 
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
        for gene in GeneSearchModel.amrList.rows 
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)

    submit: ->
        console.log("Submitting...")
        ## Checks
        if GeneSearchModel.selectedVF.length == 0 and \
           GeneSearchModel.selectedAMR.length == 0 and \
           GeneSearchModel.selectedGenomes.length == 0
            alert("You haven't selected any genes or genomes.")
        else if GeneSearchModel.selectedGenomes.length == 0
            alert("You haven't selected any genomes.")
        else if GeneSearchModel.selectedVF.length == 0 and \
                GeneSearchModel.selectedAMR.length == 0
            alert("You haven't selected any genes.")



GeneSearchModel.setLists()