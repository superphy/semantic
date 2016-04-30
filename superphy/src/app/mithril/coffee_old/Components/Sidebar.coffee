class Sidebar
	model: =>

	controller: =>

	view: =>
        m('div', {id: 'wrapper'}, [
            m('div', {id: 'sidebar-wrapper'}, [
                m('div', {class: 'sidebar'}, [
                    m('div', {id: 'sidebar-expand-collapse'}, [
                        m('button', {id: 'sidebar-expand', class: 'btn btn-default'}, [">>"])
                        m('button', {id: 'sidebar-collapse', class: 'btn btn-default'}, ["<<"])
                    ])
                ])
            ])
        ])