using EarlyBirdnBird.Data;
using Infrastructure.Interfaces;
using Infrastructure.Models;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataAccess
{
    public class UnitOfWork : IUnitOfWork
    {
        private readonly ApplicationDbContext _dbContext;

        public UnitOfWork(ApplicationDbContext dbContext)
        {
            _dbContext = dbContext;
        }

        public IGenericRepository<Amenity> _Amenity;

        public IGenericRepository<Bed> _Bed;

        public IGenericRepository<Block> _Block;

        public IGenericRepository<Image> _Image;

        public IGenericRepository<Price> _Price;

        public IGenericRepository<Rental> _Rental;

        public IGenericRepository<Reservation> _Reservation;

        public IGenericRepository<AppUser> _AppUser;

        public IGenericRepository<RentalAmenity> _RentalAmenity;

        public IGenericRepository<RentalBed> _RentalBed;

        public IGenericRepository<UserMail> _UserMail;

        public IGenericRepository<Review> _Review;



        public IGenericRepository<Amenity> Amenity
        {
            get
            {
                if (_Amenity == null)
                {
                    _Amenity = new GenericRepository<Amenity>(_dbContext);
                }
                return _Amenity;
            }
        }

        public IGenericRepository<Bed> Bed
        {
            get
            {
                if (_Bed == null)
                {
                    _Bed = new GenericRepository<Bed>(_dbContext);
                }
                return _Bed;
            }
        }

        public IGenericRepository<Block> Block
        {
            get
            {
                if (_Block == null)
                {
                    _Block = new GenericRepository<Block>(_dbContext);
                }
                return _Block;
            }
        }


        public IGenericRepository<Image> Image
        {
            get
            {
                if (_Image == null)
                {
                    _Image = new GenericRepository<Image>(_dbContext);
                }
                return _Image;
            }
        }


        public IGenericRepository<Price> Price
        {
            get
            {
                if (_Price == null)
                {
                    _Price = new GenericRepository<Price>(_dbContext);
                }
                return _Price;
            }
        }


        public IGenericRepository<Rental> Rental
        {
            get
            {
                if (_Rental == null)
                {
                    _Rental = new GenericRepository<Rental>(_dbContext);
                }
                return _Rental;
            }
        }


        public IGenericRepository<Reservation> Reservation
        {
            get
            {
                if (_Reservation == null)
                {
                    _Reservation = new GenericRepository<Reservation>(_dbContext);
                }
                return _Reservation;
            }
        }



        public IGenericRepository<AppUser> AppUser
        {
            get
            {
                if (_AppUser == null)
                {
                    _AppUser = new GenericRepository<AppUser>(_dbContext);
                }
                return _AppUser;
            }
        }


        public IGenericRepository<RentalAmenity> RentalAmenity
        {
            get
            {
                if (_RentalAmenity == null)
                {
                    _RentalAmenity = new GenericRepository<RentalAmenity>(_dbContext);
                }
                return _RentalAmenity;
            }
        }

        public IGenericRepository<RentalBed> RentalBed
        {
            get
            {
                if (_RentalBed == null)
                {
                    _RentalBed = new GenericRepository<RentalBed>(_dbContext);
                }
                return _RentalBed;
            }
        }


        public IGenericRepository<UserMail> UserMail
        {
            get
            {
                if (_UserMail == null)
                {
                    _UserMail = new GenericRepository<UserMail>(_dbContext);
                }
                return _UserMail;
            }
        }


        public IGenericRepository<Review> Review
        {
            get
            {
                if (_Review == null)
                {
                    _Review = new GenericRepository<Review>(_dbContext);
                }
                return _Review;
            }
        }


        public int Commit()
        {
            return _dbContext.SaveChanges();
        }

        public async Task<int> CommitAsync()
        {
            return await _dbContext.SaveChangesAsync();
        }

        public void Dispose()
        {
            _dbContext.Dispose();
        }
    }
}
