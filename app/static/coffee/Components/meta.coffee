#Move headers and genomes to meta class
#Follow tutorial to get asych requests 

class Singleton
    instance = null
    @get:(data)->
        instance ?= new @(data)
        return instance
        
class Meta_Data
    headers = []
    genomes = []
    meta=(data)->
        headers = data.head.vars
        for binding in data.results.bindings
            genome = []
            for variable in data.head.vars
                try
                    genome.push(binding[variable]["value"])
                catch
                    genome.push("")
            genomes.push(genome)
        {
            headers: headers
            genomes: genomes
        }

    controller: ->
        return m.request(
            method: "GET",
            url: 'http://10.139.14.121:5000/mithril/meta'
            datatype: 'json'
            type: meta)

    view: (controller) ->
        #console.log(JSON.stringify(controller))
        [
            header.view()
            m("div", {class:'container', id:'meta'},[
                m("table",{class:'well center-block'},[
                    m("tr",[
                        for item in headers
                            m("td", [item])
                    ])
                    for sequence in genomes
                        m("tr",
                            for data in sequence
                                m("td", data)
                        )
                ])
            ])
        ]