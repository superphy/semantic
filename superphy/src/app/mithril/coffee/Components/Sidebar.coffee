class Sidebar
    @controller: (args) ->
        return @
    @view: (ctrl) ->
        # m('#wrapper', [
            m('.', {id: 'sidebar-wrapper'}, [
                m('.sidebar', [
                    m('#sidebar-expand-collapse', [
                        m('button', {id: 'sidebar-expand', class: 'btn btn-default'}, [">>"])
                        m('button', {id: 'sidebar-collapse', class: 'btn btn-default'}, ["<<"])
                    ])
                ])
            ])
        # ])