## Model that holds the gene and genome data for the gene search

GeneSearchModel =
    vfs: getEndpoint(url="data/vf")
    amrs: getEndpoint(url="data/amr")
    genomes: @genomedata = getEndpoint2(url="data/meta")

    vfList: {}
    amrList: {}
    genomeList: []
    selectedVF: []
    selectedAMR: []
    ## selectedGenomes is a static list until genomes are implemented
    selectedGenomes: ["JHNV00000000", "ANVW00000000", "CP002729", "AMVC00000000"]
    ## These are setter/getters
    vfresults: m.prop({})
    amrresults: m.prop({})

    setLists: ->
        @vfList["headers"] = @vfs().headers
        @amrList["headers"] = @amrs().headers
        for gene in @vfs().rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
            @vfList["rows"] = @vfs().rows
        for gene in @amrs().rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
            @amrList["rows"] = @amrs().rows

    getSelectedVF: ->
        @selectedVF = []
        for row in @vfList.rows when row.selected()
            @selectedVF.push(row.Gene_Name)


    getSelectedAMR: ->
        @selectedAMR = []
        for row in @amrList.rows when row.selected()
            @selectedAMR.push(row.Gene_Name)

    reset: ->
        console.log("Resetting...")
        @selectedVF = []
        @selectedAMR = []
        for gene in @vfList.rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)
        for gene in @amrList.rows
            gene["selected"] = m.prop(false)
            gene["visible"] = m.prop(true)

    submit: ->
        console.log("Submitting...")
        ## Checks
        if @selectedVF.length == 0 and \
           @selectedAMR.length == 0 and \
           @selectedGenomes.length == 0
            alert("You haven't selected any genes or genomes.")
        else if @selectedGenomes.length == 0
            alert("You haven't selected any genomes.")
        else if @selectedVF.length == 0 and \
                @selectedAMR.length == 0
            alert("You haven't selected any genes.")
        else
            @vfresults = @getResults(@selectedVF)
            @amrresults = @getResults(@selectedAMR)
            m.route("/results")
            return

    getResults: (geneList) ->
        response = m.prop(null)
        if response() is null
            response =
                m.request(
                    method: "POST"
                    url: "http://#{location.hostname}:5000/data/genesearchresults",
                    data: {
                        genome: @selectedGenomes
                        genes: geneList ## temp for testing
                    }
                    datatype: "json",
                    type: (response) ->
                        console.log("Results:", response)
                        return response
                )

        return response

    getCategories: (type) ->
        response = m.prop(null)
        if response() is null
            response =
                m.request(
                    method: "POST"
                    url: "http://#{location.hostname}:5000/data/categories/#{type}"
                    data: {}
                    datatype: "json"
                    type: (response) ->
                        console.log("Categories:", response)
                        return response
                )
        return response




GeneSearchModel.setLists()