class GroupBrowse extends Page_Template
    model: () =>
        @value = m.prop('')
        @data ?= new Data()
        @table ?= new Table()
        return
    controller: (options) =>
        @model()
        return 
    view: () =>
        return [
            super()
            m('input', {oninput: m.withAttr('value', @value), value : @value() })
            m('button', { disabled : !@value(), onclick: => if @value() then @data.search(@value())},["Search"])
            @table.view(@data)
        ]