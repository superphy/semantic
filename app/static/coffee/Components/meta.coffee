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
        genomes = data.results.bindings
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
                m("table",[
                    m("tr", [
                        m("td", [item]) for item in headers
                    ])
                    for binding in genomes
                        m("tr",[ 
                            for item in headers
                                try
                                    m("td", binding[item]["value"])
                                catch
                                    m("td", "")
                        ])
                ])
            ])
        ]