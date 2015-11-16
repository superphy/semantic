urlhello = 'http://10.139.14.121:5000/mithril/query'

class Hello
    users = {}
    data = {}
    controller: () ->
        users = m.prop([])
        users = m.request({
            method: "GET", 
            url: urlhello,
            background: true, 
            initialValue: []
        }).then(users)
        users.then(m.redraw)
        
    view: () ->
        [
            home.view()
            m("div", ["Getting the first 3 triples in the triplestore."])
            #m("div",[(data.results)])
            m("p.home-beta", {class:'text-center'}, [
                JSON.stringify(users)
            ])
        ]
hello = new Hello()