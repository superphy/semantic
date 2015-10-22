gBrowse = {}

gBrowse.controller = ->

gBrowse.view = ->
    [
        header.view()
        m("div", {id: 'wrapper'}, [
            m('div', {id: 'sidebar-wrapper'})
            m('div', {id: 'page-content-wrapper'})
        ])
        sidebar.view(list.controller)
    ]
