if (window.location.hostname.includes("github.io")) {
  urlrootdir = "/OpenMS/docs"
} else if (window.location.hostname.includes("abibuilder")) {
  urlrootdir = "/archive/openms/Documentation"
} else {
  urlrootdir = "/doxygen"
}

let patharr = window.location.pathname.replace(/\/+/g, '/').split('/');

let htmlidx = patharr.indexOf("html")
let thisvers = patharr[htmlidx - 1];

$('.dropbtn').html(thisvers);
console.log(thisvers)

// https://stackoverflow.com/questions/30622369
// For GitHub Pages, use versions.json instead of directory listing
if (window.location.hostname.includes("github.io")) {
  $.getJSON(urlrootdir + '/versions.json', function(data) {
    let ret = buildVersionLinks(data.versions);
    $('.dropdown-content').append(ret.join(''));
    checkLinksExist();
  });
} else {
  $.get(urlrootdir + '/release', function(html) {
    //console.log(html)
    let ret = parseDirectoryListing(html);
    console.log("PARSED")
    $('.dropdown-content').append(ret.join(''));
    checkLinksExist();
  });
}

function checkLinksExist() {
  // Now check which links actually exist, and remove the href for those
  // that don't.
  $('.dropdown-content')
      .find('.verslink')
      .each(function () {
      var el = $(this);
      var request = new XMLHttpRequest();
      request.open('HEAD', el.attr('href'), true);
      request.onreadystatechange = function () {
          if (request.readyState === 4) {
          if (request.status === 404) {
              el.removeAttr('href');
              el.css({
              'color': "gray",
              'text-decoration': 'line-through'
              });
          }
          }
      };
      request.send();
      });
}

function buildVersionLinks(versions) {
    // Build links for GitHub Pages hosted documentation
    let docs = versions.filter(function (line) {
        // suppress link to current version
        return line != thisvers;
    });

    docs = docs.map((x) => '<a class="verslink" href="'
        + urlrootdir
        + '/' + x + '/html/' + patharr[patharr.length - 1] + '">'
        + x
        + '</a>');

    return docs;
}

function parseDirectoryListing(stringOfHtml) {
    stringOfHtml = stringOfHtml.substring(stringOfHtml.indexOf("body") - 1);
    stringOfHtml = stringOfHtml.substring(0, stringOfHtml.indexOf("/body")-1);
    var html = $($.parseHTML(stringOfHtml));
    console.log(stringOfHtml)
    var foo = []
    var docs = html.find(".fb-n>a").each(function() {
        foo.push($(this).text())
    })

    docs = foo.sort().reverse()

    docs = docs.filter(function (line) {
        // suppress link to current version
        return (line != thisvers) &&
            (line != "Parent Directory") &&
            // TODO decide how to display latest version
            (line != "latest") &&
            // suppress hidden dirs
            !/^[.]/.test(line);
    });

    docs = docs.map((x) => '<a class="verslink" href="'
        + urlrootdir
        + '/release/' + x + '/html/' + patharr[patharr.length - 1] + '">'
        + x
        + '</a>');
        
    // This assumes that a nightly documentation folder is always present.
    // A specific class does not necessarily need to be present there,
    // since this will be checked later.
    // TODO Style nightly differently to make sure users see it?
    docs.unshift('<a class="verslink" href="'
        + urlrootdir
        + '/nightly/' + 'html/' + patharr[patharr.length - 1] + '">'
        + 'nightly'
        + '</a>');
        
    return docs;
}
