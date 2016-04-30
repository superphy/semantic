# coffeelint: disable=max_line_length

class GeneData
    model: (json = {}) =>
        @categories = {}

        ## List of objects in the format of
        ## {id: , value: , title: , name: , parent: }
        @subcategories = []
        
        noID=(response)=>
            @_headers = response.head.vars
            @_rows = []
            for binding, x in response.results.bindings
                @_rows[x] = {}
                @_rows[x]["selected"] = m.prop(false)
                @_rows[x]["visible"] = m.prop(true)
                for item in response.head.vars
                    try
                        @_rows[x][item] = binding[item]['value']
                        gene_name = binding["Gene_Name"]['value']
                        category = binding["Category"]['value']
                        subcategory = binding["Sub_Category"]['value']
                        @add_gene(gene_name, category, subcategory)
                    catch
                        if item is not "Gene" then @_rows[x][item] = ''
            @search('')
        ID=(response)=>
            @_headers = ["id"].concat(response.head.vars)
            @_rows = []
            for binding, x in response.results.bindings
                @_rows[x] = {}
                @_rows[x]["id"] = x
                for item in response.head.vars
                    try
                        @_rows[x][item] = binding[item]['value']
                    catch
                        @_rows[x][item] = ''
            @search('')
        m.request(
            method: "POST",
            url: 'http://' + location.hostname + ':5000/data/' + @type,
            data: json
            datatype: 'json'
            type: noID)

    #controller
    request: (json = {}) =>
        @model(json)

    constructor: (@type, json = {}) ->
        @model(json)

    add_gene: (gene, cats, subcat) =>
        catlist = cats.split(",")
        for cat in catlist
            cat.trim()
            if cat of @categories
                if subcat of @categories[cat]
                    if gene not in @categories[cat][subcat]
                        @categories[cat][subcat].push(gene)
                else
                    @categories[cat][subcat] = [gene]
            else
                @categories[cat] = {}
                @categories[cat][subcat].push(gene)

    #view
    search: (searchterm) =>
        searchterm = searchterm.toLowerCase()
        @rows = []
        @rows.push(row) for row in @_rows when JSON.stringify(row).toLowerCase().search(searchterm) > -1
        console.log(searchterm)
        @headers = @_headers