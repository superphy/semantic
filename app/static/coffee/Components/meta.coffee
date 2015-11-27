#Move headers and genomes to meta class
#Follow tutorial to get asych requests 

class Singleton
    instance = null
    @get:(data)->
        instance ?= new @(data)
        return instance
        
class Meta_Data
    data = {}
    console.log("Meta_Data")
    meta=(response)->
        console.log("meta")
        data.headers = response.head.vars
        data.genomes = response.results.bindings
        return data

    controller: ()->
        console.log("controller")
        m.request(
            method: "GET",
            url: 'http://10.139.14.121:5000/mithril/meta'
            datatype: 'json'
            type: meta)
    views = {}
    views.table= () ->
        console.log("table")
        #console.log(JSON.stringify(controller))
        [
            m("div", {class:'container', id:'meta'},[
                m("table",[
                    m("tr", [
                        m("td", [item]) for item in data.headers
                    ])
                    for binding in data.genomes
                        m("tr",[ 
                            for item in data.headers
                                try
                                    m("td", binding[item]["value"])
                                catch
                                    m("td", "")
                        ])
                ])
            ])
        ]
    view: () ->
        console.log("view")
        [
            header.view()
            views.table()
        ]
