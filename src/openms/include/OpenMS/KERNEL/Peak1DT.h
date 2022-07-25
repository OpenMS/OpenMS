// SOA and AOS tests on god bolt

#include <tuple>
#include <vector>
#include <functional>
#include <string>
#include <limits>
#include <iostream>

template <typename T>
class proxy_holder {
  T t;
public:
  proxy_holder(T&& t) : t(t) {}
  T* operator ->() { return &t; }
};

template <typename T>
struct ref_wrap : public std::reference_wrapper<T>{
    using type = T;
    void friend swap(ref_wrap& a, ref_wrap& b) {
        T tmp = a.get();
        a.get() = b.get();
        b.get() = tmp;
     };
    void friend swap(ref_wrap&& a, ref_wrap&& b) {
        T tmp = a.get();
        a.get() = b.get();
        b.get() = tmp;
     };
    operator T&() const noexcept {return this->get(); }
    //ref_wrap(const ref_wrap& other_) : std::reference_wrapper<T>(other_.get()){}
    ref_wrap(const ref_wrap& other_) = default;
    
    //ref_wrap(ref_wrap&& other_) : std::reference_wrapper<T>(other_.get()){}
    ref_wrap(T& other_) : std::reference_wrapper<T>(other_){}
    ref_wrap(T&& other_) : std::reference_wrapper<T>(other_){}
    //void operator =(T && other_) {this->get()=other_;}
    using std::reference_wrapper<T>::operator=;
    //ref_wrap& operator=(const ref_wrap & other_) {this->get()=other_.get(); return *this;}
    //ref_wrap& operator=(ref_wrap && other_) {this->get()=other_.get(); return *this;}
    ref_wrap& operator=(T & other_) {this->get()=other_; return *this;}
    ref_wrap& operator=(T && other_) {this->get()=other_; return *this;}
};

template <typename TContainer>
class Iterator;
template <typename TContainer>
class ConstIterator;

template <template <typename...> class Container, typename TItem>
struct SOAPolicy;

template <template <typename...> class Container, template<typename...> class TItem, typename... Types>
struct SOAPolicy<Container, TItem<Types...>> {
    using type = std::tuple<Container<Types>...>;
    using const_type = std::tuple<const Container<Types>...>;
    using value_type = TItem<Types...>;
    using const_value_type = TItem<const ref_wrap<Types>...>;
    using value_ref_type = TItem<ref_wrap<Types>...>;
    using const_value_ref_type = TItem<const ref_wrap<Types>...>;
    using reference_type = TItem<ref_wrap<Types>...>&;

    constexpr static value_type get( type& c_, std::size_t position_ )
    {
        return doGet( c_, position_, std::make_integer_sequence<unsigned, sizeof...( Types )>() ); // unrolling parameter pack
    }

    constexpr static value_ref_type getRef( type& c_, std::size_t position_ )
    {
        return doGetRef( c_, position_, std::make_integer_sequence<unsigned, sizeof...( Types )>() ); // unrolling parameter pack
    }

    constexpr static const_value_ref_type getConstRef(type& c_, std::size_t position_ )
    {
        return doGetConstRef( c_, position_, std::make_integer_sequence<unsigned, sizeof...( Types )>() ); // unrolling parameter pack
    }

    constexpr static void resize( type& c_, std::size_t size_ ) {
        doResize( c_, size_, std::make_integer_sequence<unsigned, sizeof...( Types )>() ); // unrolling parameter pack
    }

    template <typename TValue>
    constexpr static void push_back( type& c_, TValue&& val_ ){
        doPushBack( c_, std::forward<TValue>(val_), std::make_integer_sequence<unsigned, sizeof...( Types )>() ); // unrolling parameter pack
    }

    static constexpr std::size_t size(type& c_){ return std::get<0>( c_ ).size(); }

    private:

    template <unsigned... Ids>
    constexpr static auto doGet( type& c_, std::size_t position_, std::integer_sequence<unsigned, Ids...> )
    {
        return value_type{ std::get<Ids>( c_ )[ position_ ]... }; // guaranteed copy elision
    }

    template <unsigned... Ids>
    constexpr static auto doGetRef( type& c_, std::size_t position_, std::integer_sequence<unsigned, Ids...> )
    {
        return value_ref_type{ ref_wrap( std::get<Ids>( c_ )[ position_ ] )... }; // guaranteed copy elision
    }

    template <unsigned... Ids>
    constexpr static auto doGetConstRef(type& c_, std::size_t position_, std::integer_sequence<unsigned, Ids...> )
    {
        return const_value_ref_type{ ref_wrap( std::get<Ids>( c_ )[ position_ ] )... }; // guaranteed copy elision
    }

    template <unsigned... Ids>
    constexpr static void doResize( type& c_, unsigned size_, std::integer_sequence<unsigned, Ids...> )
    {
        ( std::get<Ids>( c_ ).resize( size_ ), ... ); //fold expressions
    }

    template <typename TValue, unsigned... Ids>
    constexpr static void doPushBack( type& c_, TValue&& val_, std::integer_sequence<unsigned, Ids...> )
    {
        ( std::get<Ids>( c_ ).push_back( std::get<Ids>( std::forward<TValue>( val_ ) ) ), ... ); // fold expressions
    }
};


template <template <typename ValueType> class TContainer, typename TItem, typename ... ItemTs>
struct BaseContainer
{
    using policy_t        = SOAPolicy<TContainer, TItem, ItemTs ...>;
    using iterator        = Iterator<BaseContainer<TContainer, TItem>>;
    using const_iterator  = ConstIterator<BaseContainer<TContainer, TItem>>;
    using difference_type = std::ptrdiff_t;
    using value_type      = typename policy_t::value_ref_type;
    using const_value_type= typename policy_t::const_value_ref_type;
    using reference       = typename policy_t::reference_type;
    using size_type       = std::size_t;

    BaseContainer( size_t size_ ){
        resize(size_);
    }

    template <typename Fwd>
    void push_back( Fwd&& val )
    {
        policy_t::push_back( mValues, std::forward<Fwd>(val) );
    }

    std::size_t size(){
        return policy_t::size( mValues );
    }

    value_type operator[]( std::size_t position_ )
    {
        return policy_t::getRef( mValues, position_ );
    }

    const_value_type operator[]( std::size_t position_ ) const
    {
        return policy_t::getConstRef( mValues, position_ );
    }

    BaseContainer operator=(TContainer<TItem> otherAOS)
    {
        resize(0);
        for (const auto& item : otherAOS)
        {
            push_back(item);
        }
        return *this;
    }

    void resize( size_t size_ ) {
        policy_t::resize( mValues, size_ );
    }

    // Function performs insertion sort
    void insertionSort()
    {
        int N = size();
        int i, j, key;
        auto& V = *this;

        for (i = 1; i < N; i++) {
            j = i;

            // Insert V[i] into list 0..i-1
            while (j > 0 and V[j] < V[j - 1]) {

                // Swap V[j] and V[j-1]
                swap(V[j], V[j - 1]);

                // Decrement j by 1
                j -= 1;
            }
        }
    }

    iterator       begin() { return iterator( this, 0 ); }
    iterator       end() { return iterator( this, size() ); }

    // casting away const is safe here, because ConstIterator only returns const ref_wrappers
    const_iterator cbegin() const { return const_iterator( const_cast<BaseContainer*>(this), 0 ); }
    const_iterator cend() const { return const_iterator( const_cast<BaseContainer*>(this), size() ); }

    private:

    typename policy_t::type mValues;

};

template <typename TContainer>
class Iterator
{

private:
    using container_t = TContainer;
public:
    using policy_t          = typename container_t::policy_t;
    using difference_type   = std::ptrdiff_t;
    using value_type        = typename policy_t::value_ref_type;
    using reference         = typename policy_t::reference_type;
    using pointer           = value_type*;
    using iterator_category = std::random_access_iterator_tag;

    template<typename TTContainer>
    Iterator( TTContainer* container_, std::size_t position_ = 0 ):
        mContainer( container_ ),
        mIterPosition( position_ )
    {
    }

    Iterator& operator=( Iterator const& other_ ){
        mIterPosition = other_.mIterPosition;
        return *this;
        // mContainer = other_.mContainer;
    }

    difference_type operator+( Iterator const& other_ ){ std::cout << "foo+" << std::endl; return mIterPosition + other_.mIterPosition; }
    difference_type operator-( Iterator const& other_ ){ std::cout << "foo-" << std::endl; return mIterPosition - other_.mIterPosition; }
    template <typename T>
    Iterator operator+( T add ){std::cout << "foo+int" << std::endl; return Iterator(this->mContainer, this->mIterPosition + add); }
    template <typename T>
    Iterator operator-( T sub ){ return Iterator(this->mContainer, this->mIterPosition - sub); }


    friend bool operator!=( Iterator const& lhs, Iterator const& rhs ){ return lhs.mIterPosition != rhs.mIterPosition;}
    friend bool operator==( Iterator const& lhs, Iterator const& rhs ){ return !operator!=(lhs, rhs);}

    operator bool(void) const { return mIterPosition != std::numeric_limits<std::size_t>::infinity(); }

    Iterator& operator=( std::nullptr_t const& ){ mIterPosition = std::numeric_limits<std::size_t>::infinity();}

    template <typename T>
    Iterator& operator+=( T size_ ){ mIterPosition += size_; return *this;}

    template <typename T>
    Iterator& operator-=( T size_ ){ mIterPosition -= size_; return *this;}

    Iterator& operator++(){ return operator +=(1); }
    Iterator& operator--(){ return operator -=(1); }

    Iterator operator++(int){
        Iterator tmp(*this);
        operator++(); // prefix-increment this instance
        return tmp;   // return value before increment operator +=(1);
    }
    Iterator operator--(int){
        Iterator tmp(*this);
        operator--(); // prefix-decrement this instance
        return tmp;   // return value before decrement operator -=(1);
    }

    value_type operator*() {
        //std::cout << "Pos:" << mIterPosition << std::endl;
        return (*mContainer)[ mIterPosition ];
    }

    proxy_holder<value_type> operator ->() {
        return proxy_holder<value_type>(**this);
    }

    private:
    container_t*        mContainer = nullptr;
    std::size_t         mIterPosition = std::numeric_limits<std::size_t>::infinity();

};

template <typename TContainer>
class ConstIterator
{

private:
    using container_t = TContainer;
public:
    using policy_t          = typename container_t::policy_t;
    using difference_type   = std::ptrdiff_t;
    using value_type        = typename policy_t::const_value_ref_type;
    using reference         = typename policy_t::reference_type;
    using pointer           = value_type*;
    using iterator_category = std::random_access_iterator_tag;

    template<typename TTContainer>
    ConstIterator( TTContainer* container_, std::size_t position_ = 0 ):
        mContainer( container_ ),
        mIterPosition( position_ )
    {
    }

    ConstIterator& operator=( ConstIterator const& other_ ){
        mIterPosition = other_.mIterPosition;
        return *this;
        // mContainer = other_.mContainer;
    }

    difference_type operator+( ConstIterator const& other_ ){ std::cout << "foo+" << std::endl; return mIterPosition + other_.mIterPosition; }
    difference_type operator-( ConstIterator const& other_ ){ std::cout << "foo-" << std::endl; return mIterPosition - other_.mIterPosition; }
    template <typename T>
    ConstIterator operator+( T add ){std::cout << "foo+int" << std::endl; return ConstIterator(this->mContainer, this->mIterPosition + add); }
    template <typename T>
    ConstIterator operator-( T sub ){ return ConstIterator(this->mContainer, this->mIterPosition - sub); }


    friend bool operator!=( ConstIterator const& lhs, ConstIterator const& rhs ){ return lhs.mIterPosition != rhs.mIterPosition;}
    friend bool operator==( ConstIterator const& lhs, ConstIterator const& rhs ){ return !operator!=(lhs, rhs);}

    operator bool(void) const { return mIterPosition != std::numeric_limits<std::size_t>::infinity(); }

    ConstIterator& operator=( std::nullptr_t const& ){ mIterPosition = std::numeric_limits<std::size_t>::infinity();}

    template <typename T>
    ConstIterator& operator+=( T size_ ){ mIterPosition += size_; return *this;}

    template <typename T>
    ConstIterator& operator-=( T size_ ){ mIterPosition -= size_; return *this;}

    ConstIterator& operator++(){ return operator +=(1); }
    ConstIterator& operator--(){ return operator -=(1); }

    ConstIterator operator++(int){
        Iterator tmp(*this);
        operator++(); // prefix-increment this instance
        return tmp;   // return value before increment operator +=(1);
    }
    ConstIterator operator--(int){
        Iterator tmp(*this);
        operator--(); // prefix-decrement this instance
        return tmp;   // return value before decrement operator -=(1);
    }

    value_type operator*() const {
        //std::cout << "Pos:" << mIterPosition << std::endl;
        return (*mContainer)[ mIterPosition ];
    }

    proxy_holder<value_type> operator ->() const {
        return proxy_holder<value_type>(**this);
    }

    private:
    container_t*        mContainer = nullptr;
    std::size_t         mIterPosition = std::numeric_limits<std::size_t>::infinity();

};

template<typename T>
using MyVector = std::vector<T, std::allocator<T>>;

enum Component : int {
    eMz,
    eIty
};

template<typename ... T>
struct Peak1DT : public std::tuple<T ...>{
    using std::tuple<T...>::tuple;
    using std::tuple<T...>::operator=;
    Peak1DT(const Peak1DT& other) = default;
    Peak1DT(Peak1DT&& other) = default;
    Peak1DT& operator=(Peak1DT&& other_) = default;

    void friend swap(Peak1DT& a, Peak1DT& b)
    {
        a.swap(b);
    }
    void friend swap(Peak1DT&& a, Peak1DT&& b)
    {
        a.swap(b);
    }

    auto & getMZ() {return std::get<eMz>(*this);}
    auto & getMZ() const {return std::get<eMz>(*this);}
    void setMZ(double mz) {std::get<eMz>(*this) = mz;}
    auto & getIntensity() {return std::get<eIty>(*this);}
    auto & getIntensity() const {return std::get<eIty>(*this);}
    void setIntensity(double ity) {std::get<eIty>(*this) = ity;}
};

 