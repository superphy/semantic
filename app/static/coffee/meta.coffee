class Meta
	controller =
		getData: -> m.request({
			method: "GET",
			url: 'http://10.139.14.121:5000/mithril/meta'
		})
	model=
		data: m.prop(controller.getData())

	view: () ->
		[
			m("div", JSON.stringify(model.data))
		]
meta = new Meta
#class Meta
#    urlmeta = 'http://10.139.14.121:5000/mithril/meta'
#    data = {}
#    controller: () ->
#        data = m.prop([])
#        data = m.request({
#            method: "GET",
#            url: urlmeta,
#            background: true,
#            initialValue: []
#        }).then(data)
#        data.then(m.redraw)
#    view: () ->
#        [
#            home.view()
#            m("div", ["Meta-Data"])
#            #m("div",[(data.results)])
#            m("p.meta", {class:'text-center'}, [
#                JSON.stringify(JSON.parse(JSON.stringify(data)))
#            ])
#        ]
