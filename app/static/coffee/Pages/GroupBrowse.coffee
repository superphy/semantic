class GroupBrowse extends Page_Template
    model: () =>
        @value = m.prop('')
        @data ?= new Data()
        @table ?= new Table()
        return
    controller: (options) =>
        @model()
        return
    search: (searchterm) =>
        console.log(searchterm)
    view: () =>
        return [
            super()
            m('input', {oninput: m.withAttr('value', @value), value : @value() })
            m('button', { disabled : !@value(), onclick: => if @value() then @search(@value())},["Search"])
            @table.view(@data)
        ]