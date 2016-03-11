class GeneForm
    constructor: (@name, @type) ->

    model: =>

    controller: =>

    
    view: (type) =>
        m('div', {class: 'panel panel-default'}, [
            m('div', {class: 'panel-heading', id: 'vf-panel-header'}, [
                m('h4', {class: 'panel-title'}, m('a', {href: '#vf-form', config:m.route}, "#{@name} Form"))
            ])
            m('div', {class: 'panel', id: 'vf-panel'}, [
                m('div', {class: 'panel-body'}, [
                    m('div', {class: 'row'}, [
                        m('div', {class: 'col-md-6 col-md-offset-3'}, [
                            m('div', {class: 'selected-gene-list-wrapper', id: ''}'Testing Testing')
                        ])
                    ])
                ])
            ])
        ])