dragdrop = (element, options) ->
    activate = (e) ->
        e.preventDefault()
        return

    deactivate = ->

    update = (e) ->
        e.preventDefault()
        if typeof options.onchange == 'function'
            options.onchange (e.dataTransfer or e.target).files
        return

    options = options or {}
    element.addEventListener 'dragover', activate
    element.addEventListener 'dragleave', deactivate
    element.addEventListener 'dragend', deactivate
    element.addEventListener 'drop', deactivate
    element.addEventListener 'drop', update
    return


class UploadForm
    model: () =>
        @genome_name = {'value': m.prop(''), 'label':'Genome Name'}
        @serotype = {'value': m.prop(''), 'label':'Serotype'}
        @isolation_host = {'value': m.prop(''), 'label':'Isolation Host'}
        @isolation_source = {'value': m.prop(''), 'label':'Isolation Source'}

    constructor:()->
        @model()

    send: () =>
        m.request(
            method: "POST",
            url: "http://#{location.hostname}:5000/example/echo"
            data: {
                genomes: [
                    'genome_name': @genome_name.value()
                    'serotype': @serotype.value()
                    'isolation_host': @isolation_host.value()
                    'isolation_source': @isolation_source.value()
                ]
            }
            datatype: 'json',
            type: (response)=>
                console.log(response)
        )
        console.log([
            @genome_name.value()
            @serotype.value()
            @isolation_host.value()
            @isolation_source.value()
        ])
    view: () =>
        return [


            m('div', @genome_name.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @genome_name.value)
                    value: @genome_name.value()
                }))
            m('div', @serotype.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @serotype.value)
                    value: @serotype.value()
                }))
            m('div', @isolation_host.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @isolation_host.value)
                    value: @isolation_host.value()
                }))
            m('div', @isolation_source.label
                m('input[type=text]', {
                    oninput: m.withAttr('value', @isolation_source.value)
                    value: @isolation_source.value()
                }))

            m("button[type=button]", {onclick: @send}, "Done")
        ]