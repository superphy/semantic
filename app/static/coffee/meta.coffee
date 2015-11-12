urljson = 'http://10.139.14.121:5000/mithril/meta'

class Meta
	data = {}
	controller: () ->
        data = m.prop([])
        data = m.request({
            method: "GET", 
            url: urljson,
            background: true, 
            initialValue: []
        }).then(data)
        data.then(m.redraw)
        
    model : () ->
        parsed = JSON.parse(json)
        
        
    view: () ->
        [
            home.view()
            m("div", ["Meta-Data"])
            #m("div",[(data.results)])
            m("p.home-beta", {class:'text-center'}, [
                JSON.stringify(data)
            ])
        ]
meta = new Meta()