urlmeta = 'http://10.139.14.121:5000/mithril/meta'

class Meta
	data = {}
	controller: () ->
        data = m.prop([])
        data = m.request({
            method: "GET", 
            url: urlmeta,
            background: true, 
            initialValue: []
        }).then(data)
        data.then(m.redraw)
        
    view: () ->
        [
            home.view()
            m("div", ["Meta-Data"])
            #m("div",[(data.results)])
            m("p.meta", {class:'text-center'}, [
                JSON.stringify(data)
            ])
        ]
meta = new Meta()