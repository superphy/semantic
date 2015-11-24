#Move headers and genomes to meta class
#Follow tutorial to get asych requests 

class Hello
    headers = []
    genomes = []
    meta=(data)->
        setTimeout (-> say("hello world")), 50000000
        for item in data.head.vars
            headers.push(item)
        if data.results
            for i, binding in data.results.bindings
                genome = []
                for j, variable in data.head.vars
                    try
                        genome.push(data.results.bindings[binding][data.head.vars[variable]]["value"])
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
        console.log(JSON.stringify(controller))
        [
            home.view()
            m("table",{class:'text-center'},[
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
        ]